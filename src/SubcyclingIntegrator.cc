/**************************************************************************
 *
 *      Drift-mode force-subcycling integrator.
 *
 *      This is a ParaDiS port of the nonlinear-compatible mode of the
 *      ExaDiS subcycling integrator.  Segment-pair forces in groups greater
 *      than zero are cached at the current physical configuration.  Their
 *      RKF steps are virtual timestep probes and are restored afterwards.
 *      Only group zero advances the physical network, and its live force is
 *      augmented with every cached group before mobility is evaluated.
 *
 **************************************************************************/

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <vector>

#include "mpi_portability.h"

#include "Home.h"
#include "Comm.h"
#include "Subcycling.h"

namespace {

static const int MAX_SUBCYCLING_GROUPS = 5;
static const real8 MIN_SUBCYCLING_DT = 1.0e-20;

struct TagKey {
    int domain;
    int index;

    TagKey() : domain(0), index(0) {}
    explicit TagKey(const Tag_t &tag) : domain(tag.domainID), index(tag.index) {}

    bool operator<(const TagKey &other) const
    {
        return ((domain < other.domain) ||
                ((domain == other.domain) && (index < other.index)));
    }

    bool operator==(const TagKey &other) const
    {
        return ((domain == other.domain) && (index == other.index));
    }
};

struct SegmentKey {
    TagKey first;
    TagKey second;

    SegmentKey() {}

    SegmentKey(const Tag_t &tag1, const Tag_t &tag2) :
        first(tag1), second(tag2)
    {
        if (second < first) {
            TagKey tmp = first;
            first = second;
            second = tmp;
        }
    }

    bool operator<(const SegmentKey &other) const
    {
        if (first < other.first) return true;
        if (other.first < first) return false;
        return (second < other.second);
    }
};

struct SegmentPairKey {
    SegmentKey first;
    SegmentKey second;

    SegmentPairKey(Node_t *node1, Node_t *node2,
                   Node_t *node3, Node_t *node4) :
        first(node1->myTag, node2->myTag),
        second(node3->myTag, node4->myTag)
    {
        if (second < first) {
            SegmentKey tmp = first;
            first = second;
            second = tmp;
        }
    }

    bool operator<(const SegmentPairKey &other) const
    {
        if (first < other.first) return true;
        if (other.first < first) return false;
        return (second < other.second);
    }
};

struct PositionVelocity {
    Node_t *node;
    real8 x[3];
    real8 v[3];
};

struct RKFNode {
    Node_t *node;
    int native;
    real8 x[3];
    real8 v[3];
    real8 k[6][3];
};

struct RKFResult {
    real8 acceptedDT;
    real8 nextDT;
};

}  // anonymous namespace

/*
 * The definition is intentionally private to this translation unit; Home_t
 * only carries an opaque pointer while LocalSegForces uses the small helper
 * interface declared in Subcycling.h.
 */
struct _subcycling {
    int active;
    int group;
    int numGroups;
    real8 radius2[MAX_SUBCYCLING_GROUPS-1];
    real8 nextDT[MAX_SUBCYCLING_GROUPS];
    real8 realDT[MAX_SUBCYCLING_GROUPS];
    int localPairCount[MAX_SUBCYCLING_GROUPS];
    int globalPairCount[MAX_SUBCYCLING_GROUPS];
    std::map<SegmentPairKey, int> pairGroup;
    std::vector<real8> cachedForce[MAX_SUBCYCLING_GROUPS];
};

namespace {

static int TagsEqual(const Tag_t &tag1, const Tag_t &tag2)
{
    return ((tag1.domainID == tag2.domainID) &&
            (tag1.index == tag2.index));
}

static int SegmentsShareNode(Node_t *node1, Node_t *node2,
                             Node_t *node3, Node_t *node4)
{
    return (TagsEqual(node1->myTag, node3->myTag) ||
            TagsEqual(node1->myTag, node4->myTag) ||
            TagsEqual(node2->myTag, node3->myTag) ||
            TagsEqual(node2->myTag, node4->myTag));
}

static void StartSubcycling(Home_t *home)
{
    if (home->subcycling != (Subcycling_t *)NULL) {
        Fatal("Subcycling state was already active at the start of a cycle");
    }

    Param_t *param = home->param;
    Subcycling_t *state = new Subcycling_t;

    state->active = 1;
    state->group = -1;
    state->numGroups = param->subcyclingNumGroups;

    for (int i = 0; i < MAX_SUBCYCLING_GROUPS; i++) {
        state->nextDT[i] = param->subcyclingNextDT[i];
        state->realDT[i] = 0.0;
        state->localPairCount[i] = 0;
        state->globalPairCount[i] = 0;
        state->cachedForce[i].assign(3 * home->newNodeKeyPtr, 0.0);
    }

    for (int i = 0; i < state->numGroups-1; i++) {
        real8 radius = param->subcyclingRadii[i];
        state->radius2[i] = radius * radius;
    }

    /* Match the ExaDiS drift rule and explicitly catch hinges below. */
    state->radius2[0] = MAX(state->radius2[0], 1.0);
    for (int i = 1; i < state->numGroups-1; i++) {
        state->radius2[i] = MAX(state->radius2[i], state->radius2[i-1]);
    }

    home->subcycling = state;
}

static void FinishSubcycling(Home_t *home)
{
    delete home->subcycling;
    home->subcycling = (Subcycling_t *)NULL;
}

static void SetGroup(Home_t *home, int group)
{
    Subcycling_t *state = home->subcycling;
    if ((state == (Subcycling_t *)NULL) ||
        (group < 0) || (group >= state->numGroups)) {
        Fatal("Invalid subcycling force group %d", group);
    }
    state->group = group;
}

static std::vector<PositionVelocity> CapturePositionVelocity(Home_t *home)
{
    std::vector<PositionVelocity> values;
    values.reserve(home->newNodeKeyPtr + home->ghostNodeCount);

    for (int i = 0; i < home->newNodeKeyPtr; i++) {
        Node_t *node = home->nodeKeys[i];
        if (node == (Node_t *)NULL) continue;

        PositionVelocity value;
        value.node = node;
        value.x[0] = node->x;  value.x[1] = node->y;  value.x[2] = node->z;
        value.v[0] = node->vX; value.v[1] = node->vY; value.v[2] = node->vZ;
        values.push_back(value);
    }

    for (int i = 0; i < home->ghostNodeCount; i++) {
        Node_t *node = home->ghostNodeList[i];
        if (node == (Node_t *)NULL) continue;

        PositionVelocity value;
        value.node = node;
        value.x[0] = node->x;  value.x[1] = node->y;  value.x[2] = node->z;
        value.v[0] = node->vX; value.v[1] = node->vY; value.v[2] = node->vZ;
        values.push_back(value);
    }

    return values;
}

static void RestorePositionVelocity(
        const std::vector<PositionVelocity> &values)
{
    for (size_t i = 0; i < values.size(); i++) {
        Node_t *node = values[i].node;
        node->x = values[i].x[0];  node->y = values[i].x[1];
        node->z = values[i].x[2];
        node->vX = values[i].v[0]; node->vY = values[i].v[1];
        node->vZ = values[i].v[2];
    }
}

static void PreserveOuterNodalData(
        const std::vector<PositionVelocity> &values)
{
    for (size_t i = 0; i < values.size(); i++) {
        Node_t *node = values[i].node;

        node->oldx = values[i].x[0];
        node->oldy = values[i].x[1];
        node->oldz = values[i].x[2];

        node->oldvX = values[i].v[0];
        node->oldvY = values[i].v[1];
        node->oldvZ = values[i].v[2];

        node->currvX = values[i].v[0];
        node->currvY = values[i].v[1];
        node->currvZ = values[i].v[2];
    }
}

static std::vector<RKFNode> CaptureRKFNodes(Home_t *home)
{
    std::vector<RKFNode> values;
    values.reserve(home->newNodeKeyPtr + home->ghostNodeCount);

    for (int i = 0; i < home->newNodeKeyPtr; i++) {
        Node_t *node = home->nodeKeys[i];
        if (node == (Node_t *)NULL) continue;

        RKFNode value = {};
        value.node = node;
        value.native = 1;
        value.x[0] = node->x;  value.x[1] = node->y;  value.x[2] = node->z;
        value.v[0] = node->vX; value.v[1] = node->vY; value.v[2] = node->vZ;
        values.push_back(value);
    }

    for (int i = 0; i < home->ghostNodeCount; i++) {
        Node_t *node = home->ghostNodeList[i];
        if (node == (Node_t *)NULL) continue;

        RKFNode value = {};
        value.node = node;
        value.native = 0;
        value.x[0] = node->x;  value.x[1] = node->y;  value.x[2] = node->z;
        value.v[0] = node->vX; value.v[1] = node->vY; value.v[2] = node->vZ;
        values.push_back(value);
    }

    return values;
}

static void RestoreRKFBase(std::vector<RKFNode> &values)
{
    for (size_t i = 0; i < values.size(); i++) {
        Node_t *node = values[i].node;
        node->x = values[i].x[0];  node->y = values[i].x[1];
        node->z = values[i].x[2];
        node->vX = values[i].v[0]; node->vY = values[i].v[1];
        node->vZ = values[i].v[2];
    }
}

static void StoreRKFVelocity(std::vector<RKFNode> &values, int stage)
{
    for (size_t i = 0; i < values.size(); i++) {
        Node_t *node = values[i].node;
        values[i].k[stage][0] = node->vX;
        values[i].k[stage][1] = node->vY;
        values[i].k[stage][2] = node->vZ;
    }
}

static void MoveRKFNodes(Home_t *home, std::vector<RKFNode> &values,
                         real8 dt, const real8 coefficients[6])
{
    Param_t *param = home->param;

    for (size_t i = 0; i < values.size(); i++) {
        real8 x = values[i].x[0];
        real8 y = values[i].x[1];
        real8 z = values[i].x[2];

        for (int stage = 0; stage < 6; stage++) {
            x += dt * coefficients[stage] * values[i].k[stage][0];
            y += dt * coefficients[stage] * values[i].k[stage][1];
            z += dt * coefficients[stage] * values[i].k[stage][2];
        }

        FoldBox(param, &x, &y, &z);

        values[i].node->x = x;
        values[i].node->y = y;
        values[i].node->z = z;
    }
}

#ifdef ESHELBY
static void RebuildLocalIntersectionList(Home_t *home)
{
    if (home->param->enableInclusions <= 0) return;

    /*
     * Build the list for every locally owned arm at the current RKF stage.
     * This is independent of which domain happened to own a pair-force
     * calculation and therefore also makes the Eshelby mobility data robust
     * in parallel runs.
     */
    SegPartListClear(home);

    for (int i = 0; i < home->newNodeKeyPtr; i++) {
        Node_t *node = home->nodeKeys[i];
        if (node != (Node_t *)NULL) FindNodePartIntersects(home, node);
    }

    SegPartListSort(home);
}
#endif

static int GlobalMobilityError(int localError)
{
#ifdef PARALLEL
    int globalError = 0;
    MPI_Allreduce(&localError, &globalError, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    return globalError;
#else
    return localError;
#endif
}

static int EvaluateCurrentState(Home_t *home, int group, real8 dt)
{
    SetGroup(home, group);
    home->param->deltaTT = dt;

    NodeForce(home, FULL);

#ifdef ESHELBY
    /*
     * Group zero's normal base-force path has already built, exchanged,
     * and sorted the complete intersection list.  Virtual higher groups
     * skip that path, so construct the current local mobility geometry
     * explicitly for them.
     */
    if (group > 0) RebuildLocalIntersectionList(home);
#endif

    int mobilityError = CalcNodeVelocities(home, 0, 1);
    CommSendVelocity(home);

    return GlobalMobilityError(mobilityError);
}

static void SaveCurrentGroupForce(Home_t *home, int group)
{
    Subcycling_t *state = home->subcycling;
    if ((group <= 0) || (group >= state->numGroups)) return;

    std::vector<real8> &cache = state->cachedForce[group];
    size_t required = (size_t)3 * (size_t)home->newNodeKeyPtr;
    if (cache.size() != required) cache.assign(required, 0.0);

    for (int i = 0; i < home->newNodeKeyPtr; i++) {
        Node_t *node = home->nodeKeys[i];
        if (node == (Node_t *)NULL) {
            cache[3*i  ] = 0.0;
            cache[3*i+1] = 0.0;
            cache[3*i+2] = 0.0;
            continue;
        }
        cache[3*i  ] = node->fX;
        cache[3*i+1] = node->fY;
        cache[3*i+2] = node->fZ;
    }
}

static void ReduceRKFErrors(real8 localError, real8 localRelative,
                            real8 localPosition, real8 localInvalid,
                            real8 *globalError, real8 *globalRelative,
                            real8 *globalPosition, real8 *globalInvalid)
{
#ifdef PARALLEL
    real8 local[4] = { localError, localRelative,
                       localPosition, localInvalid };
    real8 global[4] = { 0.0, 0.0, 0.0, 0.0 };
    MPI_Allreduce(local, global, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    *globalError = global[0];
    *globalRelative = global[1];
    *globalPosition = global[2];
    *globalInvalid = global[3];
#else
    *globalError = localError;
    *globalRelative = localRelative;
    *globalPosition = localPosition;
    *globalInvalid = localInvalid;
#endif
}

static RKFResult AdaptiveRKFStep(Home_t *home, int group, real8 attemptDT)
{
    static const real8 stageCoefficients[6][6] = {
        { 1.0/4.0,        0.0,           0.0,            0.0,             0.0,      0.0 },
        { 3.0/32.0,       9.0/32.0,      0.0,            0.0,             0.0,      0.0 },
        { 1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0,  0.0,             0.0,      0.0 },
        { 439.0/216.0,   -8.0,           3680.0/513.0,  -845.0/4104.0,    0.0,      0.0 },
        {-8.0/27.0,       2.0,          -3544.0/2565.0,  1859.0/4104.0,  -11.0/40.0,0.0 },
        { 16.0/135.0,     0.0,           6656.0/12825.0, 28561.0/56430.0,-9.0/50.0,2.0/55.0 }
    };

    static const real8 errorCoefficients[6] = {
        1.0/360.0, 0.0, -128.0/4275.0,
       -2197.0/75240.0, 1.0/50.0, 2.0/55.0
    };

    Param_t *param = home->param;
    std::vector<RKFNode> nodes = CaptureRKFNodes(home);

    real8 dt = MIN(param->maxDT, attemptDT);
    if (dt <= 0.0) dt = param->maxDT;

    int rejected = 0;
    real8 acceptedError = 0.0;

    for (;;) {
        RestoreRKFBase(nodes);

        int mobilityError = EvaluateCurrentState(home, group, dt);

        if (!mobilityError && (group > 0)) {
            /* Cache only the force at the current physical configuration. */
            SaveCurrentGroupForce(home, group);
        }

        if (!mobilityError) {
            StoreRKFVelocity(nodes, 0);

            for (int stage = 0; stage < 5; stage++) {
                MoveRKFNodes(home, nodes, dt, stageCoefficients[stage]);
                mobilityError = EvaluateCurrentState(home, group, dt);
                if (mobilityError) break;
                StoreRKFVelocity(nodes, stage + 1);
            }
        }

        real8 globalError = 0.0;
        real8 globalRelative = 0.0;
        real8 globalPosition = 0.0;
        real8 globalInvalid = 0.0;

        if (!mobilityError) {
            real8 localError = 0.0;
            real8 localRelative = 0.0;
            real8 localPosition = 0.0;
            real8 localInvalid = 0.0;

            for (size_t i = 0; i < nodes.size(); i++) {
                if (!nodes[i].native) continue;

                real8 error[3] = { 0.0, 0.0, 0.0 };
                real8 outputDelta[3] = { 0.0, 0.0, 0.0 };

                for (int stage = 0; stage < 6; stage++) {
                    for (int component = 0; component < 3; component++) {
                        error[component] += errorCoefficients[stage] *
                                             nodes[i].k[stage][component];
                        outputDelta[component] +=
                                stageCoefficients[5][stage] *
                                nodes[i].k[stage][component];
                    }
                }

                for (int component = 0; component < 3; component++) {
                    error[component] *= dt;
                    outputDelta[component] *= dt;
                }

                real8 errorNorm = sqrt(error[0]*error[0] +
                                       error[1]*error[1] +
                                       error[2]*error[2]);
                real8 positionNorm = sqrt(outputDelta[0]*outputDelta[0] +
                                          outputDelta[1]*outputDelta[1] +
                                          outputDelta[2]*outputDelta[2]);

                localError = MAX(localError, errorNorm);
                localPosition = MAX(localPosition, positionNorm);

                real8 oldx = nodes[i].x[0];
                real8 oldy = nodes[i].x[1];
                real8 oldz = nodes[i].x[2];
                Node_t *node = nodes[i].node;
                PBCPOSITION(param, node->x, node->y, node->z,
                            &oldx, &oldy, &oldz);

                real8 dx = node->x - oldx;
                real8 dy = node->y - oldy;
                real8 dz = node->z - oldz;
                real8 displacement = sqrt(dx*dx + dy*dy + dz*dz);

                if (errorNorm > param->subcyclingRtolThreshold) {
                    real8 relative;
                    if (displacement > (param->subcyclingRtolThreshold /
                                        param->subcyclingRtolRelative)) {
                        relative = errorNorm / displacement;
                    } else {
                        relative = 2.0 * param->subcyclingRtolRelative;
                    }
                    localRelative = MAX(localRelative, relative);
                }

                if (!std::isfinite(errorNorm) ||
                    !std::isfinite(positionNorm) ||
                    !std::isfinite(displacement)) {
                    localInvalid = 1.0;
                }
            }

            ReduceRKFErrors(localError, localRelative, localPosition,
                            localInvalid, &globalError, &globalRelative,
                            &globalPosition, &globalInvalid);
        }

        int converged = (!mobilityError && (globalInvalid == 0.0) &&
                         (globalError < param->rTol) &&
                         (globalRelative < param->subcyclingRtolRelative) &&
                         (globalPosition < param->maxSeg));

        if (converged) {
            /* Fifth-order accepted position, followed by coherent forces. */
            MoveRKFNodes(home, nodes, dt, stageCoefficients[5]);
            mobilityError = EvaluateCurrentState(home, group, dt);

            if (!mobilityError) {
                acceptedError = globalError;
                break;
            }
        }

        rejected = 1;
        RestoreRKFBase(nodes);
        dt *= param->dtDecrementFact;

        if (dt < MIN_SUBCYCLING_DT) {
            Fatal("SubcyclingIntegrator(): group %d timestep dropped below "
                  "%e", group, MIN_SUBCYCLING_DT);
        }
    }

    real8 nextDT = dt;

    if (!rejected) {
        if (param->dtVariableAdjustment) {
            real8 incrementPower =
                    pow(param->dtIncrementFact, param->dtExponent);
            real8 errorRatio = acceptedError / param->rTol;
            real8 exponent = 1.0 / param->dtExponent;
            real8 factor = param->dtIncrementFact *
                    pow(1.0 / (1.0 + (incrementPower - 1.0) * errorRatio),
                        exponent);
            nextDT = dt * factor;
        } else {
            nextDT = dt * param->dtIncrementFact;
        }
    }

    RKFResult result;
    result.acceptedDT = dt;
    result.nextDT = MIN(param->maxDT, nextDT);
    return result;
}

static void ReducePairCounts(Home_t *home)
{
    Subcycling_t *state = home->subcycling;
#ifdef PARALLEL
    MPI_Allreduce(state->localPairCount, state->globalPairCount,
                  state->numGroups, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
    for (int i = 0; i < state->numGroups; i++) {
        state->globalPairCount[i] = state->localPairCount[i];
    }
#endif
}

static int SameTimeStep(real8 value1, real8 value2)
{
    real8 scale = MAX(MIN_SUBCYCLING_DT,
                      MAX(fabs(value1), fabs(value2)));
    real8 tolerance = 32.0 * std::numeric_limits<real8>::epsilon() * scale;
    return (fabs(value1 - value2) <= tolerance);
}

}  // anonymous namespace

/*------------------------------------------------------------------------
 * Public force-partition helpers
 *----------------------------------------------------------------------*/

int SubcyclingIsActive(Home_t *home)
{
    return ((home != (Home_t *)NULL) &&
            (home->subcycling != (Subcycling_t *)NULL) &&
            home->subcycling->active);
}

int SubcyclingUseBaseForces(Home_t *home)
{
    if (!SubcyclingIsActive(home)) return 1;
    return (home->subcycling->group == 0);
}

int SubcyclingSelectSegmentPair(Home_t *home, Node_t *node1, Node_t *node2,
                                Node_t *node3, Node_t *node4)
{
    if (!SubcyclingIsActive(home)) return 1;

    Subcycling_t *state = home->subcycling;
    SegmentPairKey key(node1, node2, node3, node4);
    std::map<SegmentPairKey, int>::iterator found =
            state->pairGroup.find(key);

    int group;
    if (found != state->pairGroup.end()) {
        group = found->second;
    } else {
        if (SegmentsShareNode(node1, node2, node3, node4)) {
            group = 0;
        } else {
            real8 distance2 = 0.0;
            MinSegSegDist(home, node1, node2, node3, node4, &distance2);

            group = state->numGroups - 1;
            for (int i = 0; i < state->numGroups-1; i++) {
                if (distance2 < state->radius2[i]) {
                    group = i;
                    break;
                }
            }
        }

        state->pairGroup[key] = group;
        state->localPairCount[group]++;
    }

    return (group == state->group);
}

void SubcyclingAddCachedForces(Home_t *home)
{
    if (!SubcyclingIsActive(home)) return;

    Subcycling_t *state = home->subcycling;
    if (state->group != 0) return;

    for (int i = 0; i < home->newNodeKeyPtr; i++) {
        Node_t *node = home->nodeKeys[i];
        if (node == (Node_t *)NULL) continue;

        for (int group = 1; group < state->numGroups; group++) {
            const std::vector<real8> &cache = state->cachedForce[group];
            if (cache.size() < (size_t)(3*i + 3)) continue;
            node->fX += cache[3*i  ];
            node->fY += cache[3*i+1];
            node->fZ += cache[3*i+2];
        }
    }
}

/*------------------------------------------------------------------------
 * Driver
 *----------------------------------------------------------------------*/

void SubcyclingIntegrator(Home_t *home)
{
    Param_t *param = home->param;

    StartSubcycling(home);
    Subcycling_t *state = home->subcycling;

    std::vector<PositionVelocity> outerState =
            CapturePositionVelocity(home);
    PreserveOuterNodalData(outerState);

    /*
     * The farthest group is a virtual RKF probe that selects the duration
     * of the physical outer step and initializes its cached partial force.
     */
    int highestGroup = state->numGroups - 1;
    std::vector<PositionVelocity> physicalState =
            CapturePositionVelocity(home);

    real8 highestAttempt = state->nextDT[highestGroup];
    if (highestAttempt <= 0.0) highestAttempt = param->maxDT;

    RKFResult highest = AdaptiveRKFStep(home, highestGroup,
                                        highestAttempt);
    state->realDT[highestGroup] = highest.acceptedDT;
    state->nextDT[highestGroup] = highest.nextDT;

    real8 outerDT = highest.acceptedDT;
    RestorePositionVelocity(physicalState);

    ReducePairCounts(home);

    real8 subTime[MAX_SUBCYCLING_GROUPS] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

    for (int group = 0; group < highestGroup; group++) {
        if ((group != 0) && (state->globalPairCount[group] == 0)) {
            subTime[group] = outerDT;
        }
    }

    long long numSubcycles = 0;

    while (subTime[0] < outerDT) {
        int group = highestGroup - 1;

        /* Descending scan leaves the highest-numbered group first on ties. */
        for (int candidate = highestGroup - 1;
             candidate >= 0; candidate--) {
            if (subTime[candidate] < subTime[group]) group = candidate;
        }

        real8 remaining = outerDT - subTime[group];
        if (remaining < MIN_SUBCYCLING_DT) {
            subTime[group] = outerDT;
            continue;
        }

        real8 oldNextDT = state->nextDT[group];
        if (oldNextDT <= 0.0) oldNextDT = outerDT;

        real8 attemptDT = MIN(oldNextDT, remaining);
        int clipped = (attemptDT < oldNextDT);

        std::vector<PositionVelocity> beforeVirtual;
        if (group > 0) beforeVirtual = CapturePositionVelocity(home);

        RKFResult result = AdaptiveRKFStep(home, group, attemptDT);
        state->realDT[group] = result.acceptedDT;
        state->nextDT[group] = result.nextDT;

        if (group > 0) RestorePositionVelocity(beforeVirtual);

        if (clipped && SameTimeStep(result.acceptedDT, attemptDT)) {
            /* A final clipped step must not reduce the next-cycle proposal. */
            state->nextDT[group] = oldNextDT;
        }

        subTime[group] += result.acceptedDT;
        if ((subTime[group] >= outerDT) ||
            ((outerDT - subTime[group]) < MIN_SUBCYCLING_DT)) {
            subTime[group] = outerDT;
        }

        numSubcycles++;
        if (numSubcycles > 10000000LL) {
            Fatal("SubcyclingIntegrator(): exceeded ten million subgroup "
                  "steps in one outer cycle");
        }
    }

    for (int group = 0; group < MAX_SUBCYCLING_GROUPS; group++) {
        param->subcyclingNextDT[group] = state->nextDT[group];
    }

    param->nextDT = state->nextDT[0];

    /*
     * Leave the partitioned force path before rebuilding complete final
     * nodal and arm forces.  This is required by collision/topology code.
     */
    FinishSubcycling(home);

    param->deltaTT = outerDT;
    param->realdt = outerDT;
    param->timeStart = param->timeNow;

    NodeForce(home, FULL);
    int mobilityError = CalcNodeVelocities(home, 0, 1);
    CommSendVelocity(home);

    if (GlobalMobilityError(mobilityError)) {
        Fatal("SubcyclingIntegrator(): mobility failed at the accepted final "
              "configuration");
    }

    /* CalcNodeVelocities does not alter these, but make the outer contract explicit. */
    PreserveOuterNodalData(outerState);
}
