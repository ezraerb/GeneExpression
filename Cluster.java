/* This file is part of GeneExpression, a program for plotting genes grouped by
   corrolation, the similiarity of expression behavior. It plots both a
   corrollogram and a dendrogram tree of the gene groupings

    Copyright (C) 2013   Ezra Erb

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License version 3 as published
    by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    I'd appreciate a note if you find this program useful or make
    updates. Please contact me through LinkedIn or github (my profile also has
    a link to the code depository)
*/
import java.util.*;
import java.lang.*;
import java.security.InvalidParameterException;

/* This class takes gene corrolation data and clusters genes into a dendrogram,
   a tree where each node shows the most similiar genes at that level. It
   builds the tree using the classic average-length clustering algorithm, which
   treats the corrolation values as distances between groups of genes, and
   progressively combines groups with the smallest distances. The order of
   combination defines the final dendrogram */
public final class Cluster
{
    /* Data management of the cluster distances is incredibly complex. It
       gets accessed in two very different ways: by cluster pairs during the
       combination process, and ordered by distance when finding the next pair
       to merge. Essentially, it's a sorted list which is also indexed. Note
       carefully that the two values (pair index and distance) have no
       inherent relation to each other. The sorted list will have lots of
       inserts and deletes, so this screams out for a priority queue.

       Far more accesses will be done by the cluster indexes, so the data
       structure should optimize on this operation. A multi-dimensional array
       handles it well. The priority queue then becomes a collection of array
       index pairs sorted by the array data. It gets implemented as a tree set
       with an unusual comparitor.
       
       Efficiency overview: Every merge will require updating 2N array cells,
       which requires 2N set deletes. The distance values of the new cluster
       change after the merge, so this also requires N set inserts. The set
       has N^2 entries, and both insert and delete are logorithmic, so the
       cost of either is 2logN, giving 6NlogN per merge. A dendrogram merges
       until only one cluster remains, so there will be N cluster oprations,
       which shrink by 1 with each merge, giving a final efficiency of 3NNlogN.
       Building the initial set in the first place requires a heap sort on
       N^2 entries, giving 2NNlogN, for a final cost of 5NNlogN. This is the
       same order as the theoretical minimum for average link clustering
       algorithms */

    /* Private class to hold current dendrograms representing the clusters, and
       their distances to the other clusters. Note that this is set up as a
       ragged array. Kept private to ensure consistency with the set */
    /* WARNING: Most array access operations in this object use only Java
       array checks. Its a private class, so the extra layer of security isn't
       worth it. */
    private final class ClusterDistance
    {
        private DendrogramInterface _cluster;

        /* NOTE: To avoid rounding problems, this object accululates the total
           distances for all genes in the cluster, and then calculates the
           average distance from that. The average distance is cached because
           the set operations will use it constantly */
        private double [] _totalDistances;
        private double [] _avgDistances;
        
        // Constructor. Creates a dendrogram leaf
        /* WARNING: It assumes that ALL other clusters are also dendrogram
           leaves */
        public ClusterDistance(int geneIndex, GeneCorrolations corrolations)
        {
            _cluster = new DendrogramLeaf(geneIndex, corrolations);
            /* Since the cluster contains only one gene, and all other clusters
               are leaves that contain only one gene, the distances are just
               the corrolation values for the gene */
            _totalDistances = new double[geneIndex];
            _avgDistances = new double[geneIndex];

            int index;
            for (index = 0; index < geneIndex; index++) {
                _totalDistances[index] = corrolations.getCorrolation(geneIndex,
                                                                     index);
                _avgDistances[index] = _totalDistances[index];
            }
        }

        // Get the dendrogram representing this cluster
        public DendrogramInterface getDendrogram()
        {
            return _cluster;
        }

        // Get the number of genes in this cluster
        public int getGeneCount()
        {
            return _cluster.getGeneCount();
        }

        // Extract the average distance to a given cluster
        public double getAvgDistance(int clusterIndex)
        {
            return _avgDistances[clusterIndex];
        }

        // Extract the total distance to a given cluster
        public double getTotalDistance(int clusterIndex)
        {
            return _totalDistances[clusterIndex];
        }

        // Merges another cluster into this one
        public void mergeCluster(DendrogramInterface other)
        {
            // Construct a new dendrogram branch with the two existing clusters
            DendrogramInterface newBranch = new DendrogramBranch(_cluster,
                                                                 other);
            // Overwrite the old dendrogram with the new branch
            _cluster = newBranch;
        }

        /* Accumulate the distance between this cluster and some other cluster
           due to a merge */
        public void addDistance(int clusterIndex, double otherDistance,
                                int otherSize)
        {
            _totalDistances[clusterIndex] += otherDistance;
            // To find the average distance, divide by the size of BOTH clusters
            _avgDistances[clusterIndex] = _totalDistances[clusterIndex] /
                (getGeneCount() * otherSize);
        }

        /* Clear distance data for this cluster to some other. Used when a merge
           makes a cluster disappear */
        public void clearData(int clusterIndex)
        {
            _totalDistances[clusterIndex] = 0.0;
            _avgDistances[clusterIndex] = 0.0;
        }

        // Dumps data to a string
        public String toString()
        {
            return _cluster.getGenes() + "\n Total dist:" + Arrays.toString(_totalDistances);
        }
    }

    // Corrolation data of genes to cluster
    private GeneCorrolations _corrolations;
    // Flag for whether it has been sorted and normalized
    private boolean _corrolationsSorted;
            
    // Array of cluster data
    private ClusterDistance [] _clusterData;

    // Represents a pair of clusters
    private final class ClusterPair
    {
        private int _first;
        private int _second;
        
        public ClusterPair(int first, int second)
        {
            /* Due to the ragged array, the higher index must always be
               first */
            if (first > second) {
                _first = first;
                _second = second;
            }
            else {
                _first = second;
                _second = first;
            }
        }
        
        // Getters
        public int getFirst()
        {
            return _first;
        }
        
        public int getSecond()
        {
            return _second;
        }

        // Get the distance for the cluster pair
        public double getDistance()
        {
            // NOTE: This only works due to being an inner class
            /* WARNING: If the pair is not valid, the code below can issue an
               exception. Capture it, log it, and rethrow it */
            try {
                return _clusterData[_first].getAvgDistance(_second);
            }
            catch (RuntimeException e) {
                System.out.println("Exception " + e + " caught accessing cluster pair " + toString());
                throw e;
            }
        }

        public String toString()
        {
            return "Cluster pair: " + _first + " " + _second + " dist:" +
                getDistance();
        }
    }

    // Comparitor for two cluster pairs. First compare is distance, then indexes
    private final class ClusterPairCompare implements Comparator<ClusterPair>
    {
        // Constructor. Does nothing
        public ClusterPairCompare() { ; }

        // Comparison method
        public int compare(ClusterPair pair1, ClusterPair pair2)
        {
            // Extract their distances. Lowest distance first
            double distance1 = pair1.getDistance();
            double distance2 = pair2.getDistance();

            if (distance1 < distance2)
                return -1;
            else if (distance1 > distance2)
                return 1;
            /* If distances match, sort by first cluster index (which is higher
               than the second */
            else {
                int indexDiff = pair1.getFirst() - pair2.getFirst();
                if (indexDiff != 0)
                    return indexDiff;
                else
                    // Sort by second cluster index
                    return pair1.getSecond() - pair2.getSecond();
            } // Distances match
        }

        /* Equals method for two comparitors. Since this class has no data
           members, two are equal by definition */
        public boolean equals(Object other)
        {
            return (other instanceof ClusterPairCompare);
        }
    }
        
    // Set of possible cluster pairs to merge, which manages the priority queue
    private TreeSet<ClusterPair> _possibleMergeClusters;

    // Given a pair of clusters, get the total distance
    private double getTotalClusterDistance(ClusterPair clusters)
    {
        return _clusterData[clusters.getFirst()].getTotalDistance(clusters.getSecond());
    }

    // Merges the cluster pair with the smallest average distance
    private void mergeNextClusters()
    {
        if (_possibleMergeClusters.isEmpty())
            return; // Nothing to merge!
        ClusterPair nextMerge = _possibleMergeClusters.first();

        /* The merger will invalidate all distance data for both clusters.
           Remove them all from the priority queue
           NOTE: This seems horribly inefficient. Why not just scan the
           cluster pair list? The list contains NN/2 clusters. A delete is the
           log of this, giving 2logN per delete. It's done 2N times, giving
           4NlogN. A scan requires NN comparisons, which is slower for any
           non-trivial amount of entries */
        int index;
        for (index = 0; index < _clusterData.length; index++)
            if ((_clusterData[index] != null) && // Cluster still exists
                (index != nextMerge.getFirst())) {
                _possibleMergeClusters.remove(new ClusterPair(nextMerge.getFirst(), index));
                /* NOTE: Inside the test on first to avoid trying to remove
                   nextMerge twice! */
                if (index != nextMerge.getSecond())
                    _possibleMergeClusters.remove(new ClusterPair(nextMerge.getSecond(), index));
            } // Cluster exists and not the first to merge

        // Merge the actual cluster and create the next dendrogram branch
        _clusterData[nextMerge.getFirst()].mergeCluster(_clusterData[nextMerge.getSecond()].getDendrogram());

        /* Accumulate the distances to the remaining clusters. This needs the
           size of each of those clusters as well as the distance. */
        for (index = 0; index < _clusterData.length; index++)
            if ((index != nextMerge.getFirst()) &&
                (index != nextMerge.getSecond()) &&
                (_clusterData[index] != null)) {
                ClusterPair wantPair = new ClusterPair(nextMerge.getSecond(), index);
                if (index < nextMerge.getFirst())
                    _clusterData[nextMerge.getFirst()].addDistance(index,
                                                                   getTotalClusterDistance(wantPair),
                                                                   _clusterData[index].getGeneCount());
                else
                    _clusterData[index].addDistance(nextMerge.getFirst(),
                                                    getTotalClusterDistance(wantPair),
                                                    _clusterData[index].getGeneCount());
            } // Cluster exists and not one being merged

        // Nuke the second cluster, since its now merged
        _clusterData[nextMerge.getSecond()] = null;

        // Remove data for that cluster from the ragged arrray
        for (index = nextMerge.getSecond() + 1; index < _clusterData.length;
             index++)
            if (_clusterData[index] != null)
                _clusterData[index].clearData(nextMerge.getSecond());

        /* Finally, replace the cluster indexes in the priority queue for the
           newly merged cluster. Their distances are based on the updated
           averages, so need to do last */
        for (index = 0; index < _clusterData.length; index++)
            if ((index != nextMerge.getFirst()) && (_clusterData[index] != null))
                _possibleMergeClusters.add(new ClusterPair(nextMerge.getFirst(),
                                                           index));
    }

    /* Constructor for class. This class takes ownership of the corrolations
       data and may modify it */
    public Cluster(GeneCorrolations corrolations)
    {
        _corrolations = corrolations;
        _corrolationsSorted = false;
        
        // Construct the initial cluster list, one cluster per gene
        _clusterData = new ClusterDistance[corrolations.getGeneCount()];
        int index;
        for (index = 0; index < _clusterData.length; index++)
            _clusterData[index] = new ClusterDistance(index, corrolations);
        /* Create the initial priority queue, with every possible cluster pair.
           Must be done AFTER the initial cluster data setup */

        _possibleMergeClusters = new TreeSet<ClusterPair>(new ClusterPairCompare());
        int index2;
        for (index = 0; index < _clusterData.length; index++)
            for (index2 = 0; index2 < index; index2++)
                _possibleMergeClusters.add(new ClusterPair(index, index2));
    }

    // Calculates and returns the dendrogram
    public DendrogramInterface getDendrogram()
    {
        /* NOTE: if the merge has already been done, the list of clusters to
           merge is empty and the while loop is skipped, avoiding the situation
           of calculating the dendrogram twice */
        while (!_possibleMergeClusters.isEmpty())
            mergeNextClusters();

        /* The higher cluster index gets the results after each merge, so
           the final dendrogram ends up at the end of the array */
        return _clusterData[_clusterData.length - 1].getDendrogram();
    }

    // Returns corrolation data sorted to match the dendrogram
    public GeneCorrolations getSortedCorrolations()
    {
        if (!_corrolationsSorted) {
            /* Need to sort the corrolation data based on the dendrogram. This
               requires a mapping of each gene's current position to its final
               position. Final position comes from a traversal of the
               dendrogram. Original position comes from searching the
               corrolation data by gene name. This is most efficiently done
               with a hash-map indexed by gene name */
            HashMap<String, Integer> nameMapping = new HashMap<String, Integer>();
            int index;
            for (index = 0; index < _corrolations.getGeneCount(); index++)
                nameMapping.put(_corrolations.getGeneName(index),
                                new Integer(index));

            int [] geneMapping = new int [_corrolations.getGeneCount()];

            // NOTE: The fetch call may create the dendrogram first
            int sortEntries = getSortedCorrolationsHelper(nameMapping,
                                                          geneMapping, 0,
                                                          getDendrogram());
            if (sortEntries != geneMapping.length)
                // Major, unrecoverable problem that indicates data corruption
                throw new RuntimeException();

            // Sort corrolations data
            _corrolations.reorderGenes(geneMapping);
            _corrolationsSorted = true;
        } // Sorting not already done
        return _corrolations;
    }

    /* Finds the gene sort mapping for a particular portion of a dendrogram.
       Takes the number of mappings found before it was called, and returns the
       total number done */
    private int getSortedCorrolationsHelper(HashMap<String, Integer> nameMapping,
                                            int [] geneMapping,
                                            int entriesSoFar,
                                            DendrogramInterface dendrogram)
    {
        if (dendrogram.getGeneCount() > 1) { // Branch
            entriesSoFar = getSortedCorrolationsHelper(nameMapping, geneMapping,
                                                       entriesSoFar, 
                                                       dendrogram.getFirstChild());
            entriesSoFar = getSortedCorrolationsHelper(nameMapping, geneMapping,
                                                       entriesSoFar, 
                                                       dendrogram.getSecondChild());
        }
        else { // Leaf
            /* Extract the single gene from the leaf, look it up in the gene
               name mapping, and insert that into the next entry of the gene
               index mapping */
            /* NOTE: Belive it or not, this is the cleanest method to extract
               an element from a set. If there was more than one, this would
               be a random element */
            String geneName = dendrogram.getGenes().iterator().next();
            geneMapping[entriesSoFar] = nameMapping.get(geneName).intValue();
            entriesSoFar++;
        } // Leaf
        return entriesSoFar;
    }

    // Print class data for debugging purposes
    public String toString()
    {
        StringBuffer output = new StringBuffer();
        int index;
        for (index = 0; index < _clusterData.length; index++)
            output.append(index + ": " + _clusterData[index] + '\n');
        if (!_possibleMergeClusters.isEmpty())
            output.append("Next merge: " + _possibleMergeClusters.first() + '\n');
        return output.toString();
    }

    // Object test program. Reads a hard coded data file and clusters it
    public static void main(String[] args) throws Throwable
    {
        try {
            DataLoader data = new DataLoader();
            data.loadData("data", "test.txt", 5.0);
            GeneCorrolations geneData = new GeneCorrolations(data);
            Cluster test = new Cluster(geneData);
            System.out.println(test);
            
            // Cluster the data, tracing state after each step
            while (!test._possibleMergeClusters.isEmpty()) {
                test.mergeNextClusters();
                System.out.println(test);
            }

            // Print the final dendrogram
            System.out.println(test.getDendrogram());
            System.out.println(test.getSortedCorrolations());
        }
        catch (Throwable e) {
            System.out.println("Exeption " + e + " caught");
            throw e; // Force improper termination, so error is obvious
        }
    }
}
