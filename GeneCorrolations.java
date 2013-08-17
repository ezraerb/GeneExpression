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

/* This class calculates and holds similarity data of gene expression patterns.
   It uses the Ecludian distance, the most widely used metric on normalized
   data */
public final class GeneCorrolations
{
    /* Gene corrolation is a map of every gene to every other gene. To minimize
       space usage, it is held as a ragged array. This class has a method to
       access the array as an actual two dimensional array, but it is less
       efficient than scanning the array directly */

    // Private class to hold gene corrolations
    private final class Corrolate
    {
        private String _geneName;
        private double [] _corrolations;

        // Constructor
        public Corrolate(String geneName, double [] corrolations)
        {
            _geneName = geneName;
            _corrolations = corrolations;
        }

        // Getters
        public String getGeneName()
        {
            return _geneName;
        }

        // Extract the corolation to a given gene
        public double getCorrolate(int gene)
        {
            if (gene >= _corrolations.length)
                throw new InvalidParameterException();
            return _corrolations[gene];
        }

        // Returns the maximum corrolation value for this set of data
        /* NOTE: This is not necessarily the highest value for the gene, thanks
           to the ragged array setup */
        public double maxValue()
        {
            /* Believe it or not, Java has no standard function to find the
               maximum value in an array, leading to this design pattern */
            double maxValue = 0.0;
            int index;
            for (index = 0; index < _corrolations.length; index++)
                if (_corrolations[index] > maxValue)
                    maxValue = _corrolations[index];
            return maxValue;
        }
        
        // Normalize the corrolation values given the max corrolation value
        public void normalize(double maxValue)
        {
            int index;
            for (index = 0; index < _corrolations.length; index++)
                _corrolations[index] /= maxValue;
        }

        // Dumps data to a string
        public String toString()
        {
            return _geneName + ":" + Arrays.toString(_corrolations);
        }
    }

    private Corrolate[] _corrolations;

    /* Constructor. Builds the corrolation matrix based on the passed in
       file data */
    public GeneCorrolations(DataLoader data)
    {
        // Get a pointer to the raw gene data, to avoid extracting it constantly
        ArrayList<GeneData> geneData = data.getGeneData();

        _corrolations = new Corrolate[geneData.size()];

        /* The corrolation matrix requires computing every gene with every other
           gene, requiring a double loop. Since the array is ragged, the row
           index also gives the size of the row */

        double [] results = null;
        int index1, index2;
        for (index1 = 0; index1 < _corrolations.length; index1++) {
            results = new double[index1];
            for (index2 = 0; index2 < index1; index2++)
                results[index2] = findCorrolation(geneData.get(index1).getGeneData(),
                                                  geneData.get(index2).getGeneData());
            _corrolations[index1] = new Corrolate(geneData.get(index1).getGeneName(),
                                                  results);
        } // Outer loop over genes

        // Normalize the data
        normalize();
    } // Constructor

    // Finds the corrolation between two genes given observation data
    private double findCorrolation(double [] first, double [] second)
    {
        /* The ecludian distane is defined as the distance between the two
           data points in an N-dimensional space */
        if ((first.length != second.length) || (first.length == 0))
            // Inconsistent data
            throw new InvalidParameterException();
        double result = 0.0; 
        int index;
        for (index = 0; index < first.length; index++)
            result += (first[index] - second[index]) * (first[index] - second[index]);
        result /= first.length;
        return Math.sqrt(result);
    }

    // Normalizes corrolations to be in the range 0..1
    private void normalize()
    {
        /* To normalize the corrolation values, first need to find the highest
           corrolation value. Its expensive but straightforward */
        double maxValue = 0.0;
        int index;
        for (index = 0; index < _corrolations.length; index++) {
            double newMax = _corrolations[index].maxValue();
            if (newMax > maxValue)
                maxValue = newMax;
        }
        // If the value is greater than zero, divide all values by this value
        if (maxValue > 0.0)
            for (index = 0; index < _corrolations.length; index++)
                _corrolations[index].normalize(maxValue);
    }

    // Reorders the corrolation matrix by the specified mapping
    public void reorderGenes(int [] newOrdering)
    {
        /* The reordering list shows for position, the index of the existing
           matrix to move to that spot. Need to confirm all of the following:
           1. It has the right number of entries.
           2. Each of its values is valid
           3. None are duplicated
           The test is not cheap, but worth it because an invalid ordering
           will invalidate the matrix in ways that are hard to detect */
        if (_corrolations.length != newOrdering.length)
            throw new InvalidParameterException();
        /* To check for duplicates, create a set of the values. If the final
           size does not equal the array size, have a duplicate somewhere */
        Set<Integer> allValues = new HashSet<Integer>();

        int index;
        for (index = 0; index < newOrdering.length; index++) {
            if ((newOrdering[index] < 0) ||
                (newOrdering[index] >= _corrolations.length))
                throw new InvalidParameterException();
            else
                allValues.add(new Integer(newOrdering[index]));
        } // Loop through the map values
        if (allValues.size() != newOrdering.length)
            throw new InvalidParameterException();

        // Create a new corrolation matrix, and then replace the existing one
        Corrolate [] newMatrix = new Corrolate[_corrolations.length];
        double [] newData = null;

        int index1, index2;
        for (index1 = 0; index1 < newOrdering.length; index1++) {
            newData = new double[index1];
            for (index2 = 0; index2 < index1; index2++)
                newData[index2] = getCorrolation(newOrdering[index1],
                                                 newOrdering[index2]);
            newMatrix[index1] = new Corrolate(getGeneName(newOrdering[index1]),
                                              newData);
        } // For loop on corrolate matrix rows
        _corrolations = newMatrix;
    }

    // Fetch the gene name given it's index in the matrix
    public String getGeneName(int index)
    {
        if ((index < 0) || (index >= _corrolations.length))
            throw new InvalidParameterException();
        return _corrolations[index].getGeneName();
    }

    // Fetch a corrolation value given the indexes of the two genes
    public double getCorrolation(int gene1, int gene2)
    {
        if ((gene1 < 0) || (gene1 >= _corrolations.length) ||
            (gene2 < 0) || (gene2 >= _corrolations.length))
            throw new InvalidParameterException();
        // Due to the ragged array, look up the higher index first
        else if (gene1 > gene2)
            return _corrolations[gene1].getCorrolate(gene2);
        else if (gene2 > gene1)
            return _corrolations[gene2].getCorrolate(gene1);
        else // Corrolation with self is always zero
            return 0.0;
    }

    // Returns the number of genes with corrolations defined
    public int getGeneCount()
    {
        return _corrolations.length;
    }

    // Dumps the contects of the class to a string
    public String toString()
    {
        StringBuffer output = new StringBuffer();
        int index;
        for (index = 0; index < _corrolations.length; index++) {
            output.append(_corrolations[index]);
            output.append('\n');
        }
        return output.toString();
    }

    // Object test program. Reads a hard coded data file, outputs the results
    public static void main(String[] args) throws Throwable
    {
        try {
            DataLoader data = new DataLoader();
            data.loadData("data", "test.txt", 5.0);
            GeneCorrolations test = new GeneCorrolations(data);
            int [] testMapping = new int[test.getGeneCount()];

            // Flip the gene order
            int index;
            for (index = 0; index < testMapping.length; index++)
                testMapping[index] = testMapping.length - index - 1;
            test.reorderGenes(testMapping);

            System.out.println(test);
        }
        catch (Throwable e) {
            System.out.println("Exeption " + e + " caught");
            throw e; // Force improper termination, so error is obvious
        }
    }
}