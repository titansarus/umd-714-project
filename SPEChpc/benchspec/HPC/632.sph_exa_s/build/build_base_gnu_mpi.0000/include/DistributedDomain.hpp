#pragma once

#ifdef USE_MPI
#include "mpi.h"
#endif

#include <vector>
#include <cmath>
#include <map>
#include <algorithm>

#include "BBox.hpp"

namespace sphexa
{

template <typename T>
class DistributedDomain
{
#ifdef USE_MPI
public:
    DistributedDomain(MPI_Comm comm = MPI_COMM_WORLD)
        : comm(comm)
    {
        MPI_Comm_size(comm, &comm_size);
        MPI_Comm_rank(comm, &comm_rank);
        MPI_Get_processor_name(processor_name, &name_len);
    }

    inline T normalize(T d, T min, T max) { return (d - min) / (max - min); }

    inline bool overlap(T leftA, T rightA, T leftB, T rightB) { return leftA < rightB && rightA > leftB; }

    inline T computeMaxH(const std::vector<int> &clist, const Array<T> &h)
    {
        T hmax = 0.0;
        for (unsigned int i = 0; i < clist.size(); i++)
        {
            T hh = h[clist[i]];
            if (hh > hmax) hmax = hh;
        }
        return hmax;
    }

    inline T computeGlobalMaxH(const std::vector<int> &clist, const Array<T> &h)
    {
        T hmax = computeMaxH(clist, h);

        MPI_Allreduce(MPI_IN_PLACE, &hmax, 1, MPI_DOUBLE, MPI_MAX, comm);

        return hmax;
    }

    void distributeParticles(const std::vector<int> &clist, const BBox<T> &globalBBox, const Array<T> &x, const Array<T> &y,
                             const Array<T> &z)
    {
        cellList.clear();
        cellList.resize(ncells);

        std::vector<int> localCount(ncells);

        for (unsigned int i = 0; i < clist.size(); i++)
        {
            T xx = std::max(std::min(x[clist[i]], globalBBox.xmax), globalBBox.xmin);
            T yy = std::max(std::min(y[clist[i]], globalBBox.ymax), globalBBox.ymin);
            T zz = std::max(std::min(z[clist[i]], globalBBox.zmax), globalBBox.zmin);

            T posx = normalize(xx, globalBBox.xmin, globalBBox.xmax);
            T posy = normalize(yy, globalBBox.ymin, globalBBox.ymax);
            T posz = normalize(zz, globalBBox.zmin, globalBBox.zmax);

            int hx = posx * nX;
            int hy = posy * nY;
            int hz = posz * nZ;

            hx = std::min(hx, nX - 1);
            hy = std::min(hy, nY - 1);
            hz = std::min(hz, nZ - 1);

            unsigned int l = hz * nX * nY + hy * nX + hx;

            cellList[l].push_back(clist[i]);
        }
    }

    void computeGlobalCellCount()
    {
        globalCellCount.resize(ncells);

        std::vector<int> localCount(ncells);

        for (int i = 0; i < ncells; i++)
            localCount[i] = cellList[i].size();

        MPI_Allreduce(&localCount[0], &globalCellCount[0], ncells, MPI_INT, MPI_SUM, comm);
    }

    void assignRanks(const std::vector<int> &procsize)
    {
        assignedRanks.resize(ncells);
        for (unsigned int i = 0; i < assignedRanks.size(); i++)
            assignedRanks[i] = 0;

        int rank = 0;
        int work = 0;

        // Assign rank to each cell
        for (int i = 0; i < ncells; i++)
        {
            work += globalCellCount[i];
            assignedRanks[i] = rank;

            if (work > 0 && work >= procsize[rank])
            {
                work = 0;
                rank++;
            }
        }
    }

    inline void resize(unsigned int size, std::vector<Array<T> *> &data)
    {
        for (unsigned int i = 0; i < data.size(); i++)
            data[i]->resize(size);
    }

    void synchronize(std::vector<Array<T> *> &data)
    {
        std::map<int, std::vector<int>> toSend;
        std::vector<std::vector<T>> buff;
        std::vector<int> counts;

        int needed = 0;

        for (int i = 0; i < ncells; i++)
        {
            int rank = assignedRanks[i];
            if (rank != comm_rank && cellList[i].size() > 0)
                toSend[rank].insert(toSend[rank].end(), cellList[i].begin(), cellList[i].end());
            else if (rank == comm_rank)
                needed += globalCellCount[i] - cellList[i].size();
        }

        std::vector<MPI_Request> requests;
        counts.resize(comm_size);
        for (int rank = 0; rank < comm_size; rank++)
        {
            if (toSend[rank].size() > 0)
            {
                int rcount = requests.size();
                int bcount = buff.size();
                counts[rank] = toSend[rank].size();

                requests.resize(rcount + data.size() + 2);
                buff.resize(bcount + data.size());

                MPI_Isend(&comm_rank, 1, MPI_INT, rank, 0, comm, &requests[rcount]);
                MPI_Isend(&counts[rank], 1, MPI_INT, rank, 1, comm, &requests[rcount + 1]);

                for (unsigned int i = 0; i < data.size(); i++)
                {
                    buff[bcount + i].resize(counts[rank]);

                    for (int j = 0; j < counts[rank]; j++)
                        buff[bcount + i][j] = (*data[i])[toSend[rank][j]];

                    MPI_Isend(&buff[bcount + i][0], counts[rank], MPI_DOUBLE, rank, 2 + i, comm, &requests[rcount + 2 + i]);
                }
            }
        }

        int end = data[0]->size();
        resize(end + needed, data);

        while (needed > 0)
        {
            std::vector<MPI_Status> status(2);

            int rank, count;
            MPI_Recv(&rank, 1, MPI_INT, MPI_ANY_SOURCE, 0, comm, &status[0]);
            MPI_Recv(&count, 1, MPI_INT, rank, 1, comm, &status[1]);

            int rcount = requests.size();
            requests.resize(rcount + data.size());

            for (unsigned int i = 0; i < data.size(); i++)
                MPI_Irecv(&(*data[i])[end], count, MPI_DOUBLE, rank, 2 + i, comm, &requests[rcount + i]);


            end += count;
            needed -= count;
        }

        if (requests.size() > 0)
        {
            std::vector<MPI_Status> status(requests.size());
            MPI_Waitall(requests.size(), &requests[0], &status[0]);
        }
    }

    bool checkBoxOverlap(const std::vector<BBox<T>> &cellBBox, T xmin, T xmax, T ymin, T ymax, T zmin, T zmax, T ri, int l)
    {
        return overlap(cellBBox[l].xmin, cellBBox[l].xmax, xmin - ri, xmax + ri) &&
               overlap(cellBBox[l].ymin, cellBBox[l].ymax, ymin - ri, ymax + ri) &&
               overlap(cellBBox[l].zmin, cellBBox[l].zmax, zmin - ri, zmax + ri);
    }

    void computeHaloList(const BBox<T> &localBBox, const BBox<T> &globalBBox, bool showGraph = false)
    {
        sendHaloList.clear();
        recvHaloList.clear();

        int reorder = 0;
        haloCount = 0;
        MPI_Comm graphComm;
        std::map<int, std::vector<int>> msources;

        {
            std::vector<bool> visited(ncells, false);
            std::vector<int> sources, weights, degrees, dests;

            int mix = (int)floor(normalize(localBBox.xmin - 2 * localMaxH, globalBBox.xmin, globalBBox.xmax) * nX);
            int miy = (int)floor(normalize(localBBox.ymin - 2 * localMaxH, globalBBox.ymin, globalBBox.ymax) * nY);
            int miz = (int)floor(normalize(localBBox.zmin - 2 * localMaxH, globalBBox.zmin, globalBBox.zmax) * nZ);
            int max = (int)floor(normalize(localBBox.xmax + 2 * localMaxH, globalBBox.xmin, globalBBox.xmax) * nX);
            int may = (int)floor(normalize(localBBox.ymax + 2 * localMaxH, globalBBox.ymin, globalBBox.ymax) * nY);
            int maz = (int)floor(normalize(localBBox.zmax + 2 * localMaxH, globalBBox.zmin, globalBBox.zmax) * nZ);

            if (!globalBBox.PBCx) mix = std::max(mix, 0);
            if (!globalBBox.PBCy) miy = std::max(miy, 0);
            if (!globalBBox.PBCz) miz = std::max(miz, 0);
            if (!globalBBox.PBCx) max = std::min(max, nX - 1);
            if (!globalBBox.PBCy) may = std::min(may, nY - 1);
            if (!globalBBox.PBCz) maz = std::min(maz, nZ - 1);

            for (int hz = miz; hz <= maz; hz++)
            {
                for (int hy = miy; hy <= may; hy++)
                {
                    for (int hx = mix; hx <= max; hx++)
                    {
                        // T displz = bbox.PBCz? ((hz < 0) - (hz >= nZ)) * (globalBBox.zmax-globalBBox.zmin) : 0;
                        // T disply = bbox.PBCy? ((hy < 0) - (hy >= nY)) * (globalBBox.ymax-globalBBox.ymin) : 0;
                        // T displx = bbox.PBCx? ((hx < 0) - (hx >= nX)) * (globalBBox.xmax-globalBBox.xmin) : 0;

                        int hzz = globalBBox.PBCz ? (hz % nZ) + (hz < 0) * nZ : hz;
                        int hyy = globalBBox.PBCy ? (hy % nY) + (hy < 0) * nY : hy;
                        int hxx = globalBBox.PBCx ? (hx % nX) + (hx < 0) * nX : hx;

                        unsigned int l = hzz * nY * nX + hyy * nX + hxx;

                        if (!visited[l] && assignedRanks[l] != comm_rank &&
                            globalCellCount[l] >
                                0) // && checkBoxOverlap(cellBBox, localBBox.xmin+displx, localBBox.xmax+displx, localBBox.ymin+disply,
                                   // localBBox.ymax+disply, localBBox.zmin+displz, localBBox.zmax+displz, 2*localMaxH, l))
                        {
                            int rank = assignedRanks[l];
                            msources[rank].push_back(l);
                            recvHaloList[rank] += globalCellCount[l];
                            haloCount += globalCellCount[l];
                            visited[l] = true;
                        }
                    }
                }
            }

            int n = msources.size();
            sources.resize(n);
            weights.resize(n);
            degrees.resize(n);
            dests.resize(n);

            int i = 0;
            for (auto it = msources.begin(); it != msources.end(); ++it)
            {
                sources[i] = it->first;
                weights[i] = it->second.size();
                degrees[i] = 1;
                dests[i++] = comm_rank;
            }

            MPI_Dist_graph_create(comm, n, &sources[0], &degrees[0], &dests[0], &weights[0], MPI_INFO_NULL, reorder, &graphComm);
        }

        int indegree = 0, outdegree = 0, weighted = 0;
        MPI_Dist_graph_neighbors_count(graphComm, &indegree, &outdegree, &weighted);

        std::vector<int> sources(indegree), sourceweights(indegree);
        std::vector<int> dests(outdegree), destweights(outdegree);

        MPI_Dist_graph_neighbors(graphComm, indegree, &sources[0], &sourceweights[0], outdegree, &dests[0], &destweights[0]);

        if (showGraph)
        {
            printf("\t%d -> { ", comm_rank);
            for (int i = 0; i < indegree; i++)
                printf("%d ", sources[i]);
            printf("};\n");
            fflush(stdout);
        }

        std::vector<MPI_Request> requests(indegree);

        // Ask sources a list of cells
        // msources contains a list of cell id for each rank
        for (int i = 0; i < indegree; i++)
        {
            int rank = sources[i];
            int count = msources[rank].size();
            // printf("%d send %d to %d\n", comm_rank, count, rank); fflush(stdout);
            MPI_Isend(&msources[rank][0], count, MPI_INT, rank, 0, graphComm, &requests[i]);
        }

        // Recv cell list from destinations
        // buff contains a list of cell that we will have to send
        // senHaloList contains for each rank, the list of particles to cell
        // Better to move this to a list of cells
        for (int i = 0; i < outdegree; i++)
        {
            MPI_Status status;
            int rank = dests[i];
            int count = destweights[i];
            std::vector<int> buff(count);
            // printf("%d rcv %d from %d\n", comm_rank, count, rank); fflush(stdout);
            MPI_Recv(&buff[0], count, MPI_INT, rank, 0, graphComm, &status);
            for (int j = 0; j < count; j++)
                sendHaloList[rank].insert(sendHaloList[rank].end(), cellList[buff[j]].begin(), cellList[buff[j]].end());
        }

        if (requests.size() > 0)
        {
            std::vector<MPI_Status> status(requests.size());
            MPI_Waitall(requests.size(), &requests[0], &status[0]);
        }

        MPI_Comm_free(&graphComm);
    }

    void makeDataArray(std::vector<Array<T> *> &data, Array<T> *d) { data.push_back(d); }

    template <typename... Args>
    void makeDataArray(std::vector<Array<T> *> &data, Array<T> *first, Args... args)
    {
        data.push_back(first);
        makeDataArray(data, args...);
    }

    template <typename... Args>
    void synchronizeHalos(Args... args)
    {
        std::vector<Array<T> *> data;
        makeDataArray(data, args...);
        synchronizeHalos(data);
    }

    void synchronizeHalos(std::vector<Array<T> *> &data)
    {
        std::vector<std::vector<T>> buff;
        std::vector<MPI_Request> requests;

        for (auto it = sendHaloList.begin(); it != sendHaloList.end(); ++it)
        {
            int rank = it->first;
            const std::vector<int> &cellList = it->second;

            int rcount = requests.size();
            int bcount = buff.size();
            int count = cellList.size();

            requests.resize(rcount + data.size());
            buff.resize(bcount + data.size());

            for (unsigned int i = 0; i < data.size(); i++)
            {
                buff[bcount + i].resize(count);

                for (int j = 0; j < count; j++)
                    buff[bcount + i][j] = (*data[i])[cellList[j]];

                MPI_Isend(&buff[bcount + i][0], count, MPI_DOUBLE, rank, i, comm, &requests[rcount + i]);
            }
        }

        int end = data[0]->size();
        resize(end + haloCount, data);

        for (auto it = recvHaloList.begin(); it != recvHaloList.end(); ++it)
        {
            std::vector<MPI_Status> status(data.size());

            int rank = it->first;
            int count = it->second;

            for (unsigned int i = 0; i < data.size(); i++)
                MPI_Recv(&(*data[i])[end], count, MPI_DOUBLE, rank, i, comm, &status[i]);

            end += count;
        }

        if (requests.size() > 0)
        {
            std::vector<MPI_Status> status(requests.size());
            MPI_Waitall(requests.size(), &requests[0], &status[0]);
        }
    }

    inline void removeIndices(const std::vector<bool> indices, std::vector<Array<T> *> &data)
    {
        for (unsigned int i = 0; i < data.size(); i++)
        {
            Array<T> &array = *data[i];

            int j = 0;
            std::vector<T> tmp(array.size());
            for (unsigned int i = 0; i < array.size(); i++)
            {
                if (indices[i] == false) tmp[j++] = array[i];
            }
            tmp.swap(array);
            array.resize(j);
        }
    }

    void computeDiscardList(const int count, std::vector<bool> &discardList)
    {
        discardList.resize(count, false);
        for (int i = 0; i < ncells; i++)
        {
            if (assignedRanks[i] != comm_rank)
            {
                for (unsigned j = 0; j < cellList[i].size(); j++)
                    discardList[cellList[i][j]] = true;
            }
        }
    }

    void discard(std::vector<Array<T> *> &data)
    {
        std::vector<bool> discardList;
        computeDiscardList(data[0]->size(), discardList);
        if (discardList.size() > 0) removeIndices(discardList, data);
    }

    template <class Dataset>
    void distribute(std::vector<int> &clist, Dataset &d, bool showGraph = false)
    {
        /* The 'bbox' here is is only used to test PBC. If PBC is activated, the globalBBox will take the corresponding bbox.x{min,max}
         * values */
        d.bbox.computeGlobal(clist, d.x, d.y, d.z, comm);
        globalMaxH = computeGlobalMaxH(clist, d.h);

        nX = std::max((d.bbox.xmax - d.bbox.xmin) / globalMaxH, 2.0);
        nY = std::max((d.bbox.ymax - d.bbox.ymin) / globalMaxH, 2.0);
        nZ = std::max((d.bbox.zmax - d.bbox.zmin) / globalMaxH, 2.0);
        ncells = nX * nY * nZ;

        distributeParticles(clist, d.bbox, d.x, d.y, d.z); // Distribute particles in global buckets
        computeGlobalCellCount();                          // How many particles per bucket (globally) ?
        assignRanks(d.workload);                           // Assign ranks to buckets
        synchronize(d.data);                               // Receive new particles from neighbor nodes
        discard(d.data);                                   // Discard particles that were sent to other nodes

        clist.resize(d.data[0]->size());
        for (unsigned int i = 0; i < clist.size(); i++)
            clist[i] = i;

        distributeParticles(clist, d.bbox, d.x, d.y, d.z); // Re-distribute new particles in global buckets

        BBox<T> localBBox;
        localBBox.compute(clist, d.x, d.y, d.z);
        localMaxH = computeMaxH(clist, d.h);

        // Use the localbbox to identify halo cells
        computeHaloList(localBBox, d.bbox, showGraph); // Find which rank is our neighbor (builds a directed graph)

        d.count = clist.size(); // Update particle count
    }

private:
    MPI_Comm comm;

    T localMaxH, globalMaxH;

    std::vector<int> assignedRanks;
    std::vector<int> globalCellCount;
    std::vector<std::vector<int>> cellList;

    std::map<int, std::vector<int>> sendHaloList;
    std::map<int, int> recvHaloList;

    int ncells;
    int nX, nY, nZ;

    int comm_size, comm_rank, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

#else
public:
    template <class Dataset>
    void distribute(const std::vector<int> &clist, Dataset &d)
    {
        d.bbox.compute(clist, d.x, d.y, d.z);
    }

    template <typename... Args>
    void synchronizeHalos(Args...)
    {
    }
#endif

public:
    template <typename... Args>
    void resizeArrays(const int count, Array<T> *d)
    {
        d->resize(count);
    }

    template <typename... Args>
    void resizeArrays(const int count, Array<T> *first, Args... args)
    {
        first->resize(count);
        resizeArrays(count, args...);
    }

public:
    int haloCount = 0;
};

} // namespace sphexa
