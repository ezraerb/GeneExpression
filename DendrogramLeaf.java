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
import java.security.InvalidParameterException;

/* This class defines the leaf of a dendrogram, a tree that groups genes by
   similarity. Its implemented using the classic whole-part design pattern,
   so all nodes act identically */
public final class DendrogramLeaf implements DendrogramInterface
{
    /* The single gene in the leaf. Stored as a set to match the interface
       without calling set constructors everywhere */
    private Set<String> _gene;
    
    /* Node weight, defined as the total of all expression similiarities for
       this gene */
    double _weight;
    
    /* Constructor. Needs the index of the gene in the corrolations map,
       and the corrolation data */
    public DendrogramLeaf(int geneIndex, GeneCorrolations corrolations)
    {
        if (corrolations == null)
            throw new InvalidParameterException();
        if ((geneIndex < 0) || (geneIndex >= corrolations.getGeneCount()))
            throw new InvalidParameterException();
        _gene = new HashSet<String>();
        _gene.add(corrolations.getGeneName(geneIndex));
        
        /* Get the weight, defined as the sum of similiarities to other genes.
           Lower weight is better */
        /* TRICKY NOTE: Due to the ragged array within the corrolations data,
           calling with the higher gene index first will be faster */
        _weight = 0.0;
        int index;
        for (index = 0; index < geneIndex; index++)
            _weight += corrolations.getCorrolation(geneIndex, index);
        for (index = (geneIndex + 1); index < corrolations.getGeneCount();
             index++)
            _weight += corrolations.getCorrolation(index, geneIndex);
    }

    /* Number of genes in this node or its children. A value of one indicates
       a leaf */
    public int getGeneCount()
    {
        return 1; // By definition
    }

    // Name of gene in the node
    public Set<String> getGenes()
    {
        return _gene;
    }

    /* Total of the weigts of all leaves under this node, or the weight if
       its a leaf */
    public double getTotalWeight()
    {
        return _weight;
    }

    /* Level of this node, defined as the number of nodes to reach the
       furthest leaf */
    public int getLevel()
    {
        return 0; // True by definition
    }
    
    // Returns first child, if any
    public DendrogramInterface getFirstChild()
    {
        return null; // Leaves have no children by definition
    }
    
    // Returns second child, if any
    public DendrogramInterface getSecondChild()
    {
        return null;
    }

    /* Averegate weight of the node, defined as the average leaf weight, used
       for sorting */
    public double getAverageWeight()
    {
        // Every leaf holds a single gene, so total weight is also the average
        return getTotalWeight();
    }

    // Output the portion of the tree of this node and below
    public void toString(StringBuffer output, int maxLevel, int parentLevel)
    {
        int index;
        if (maxLevel > getLevel()) {
            if (parentLevel > getLevel()) { // This node has a parent
                for (index = maxLevel; index > parentLevel; index--)
                    output.append(' ');
                for (index = parentLevel; index > getLevel(); index--)
                    output.append('-');
            }
            else
                for (index = maxLevel; index > getLevel(); index--)
                    output.append(' ');
        }
        output.append("LEAF Gene: " + _gene + " weight: " + _weight + '\n');
    }
    
    // Output the tree starting with this node
    public String toString()
    {
        StringBuffer output = new StringBuffer();
        toString(output, getLevel(), 0);
        return output.toString();
    }

   // Object test program. Reads a hard coded data file, outputs the results
    public static void main(String[] args) throws Throwable
    {
        try {
            DataLoader data = new DataLoader();
            data.loadData("data", "test.txt", 5.0);
            GeneCorrolations geneData = new GeneCorrolations(data);
            DendrogramLeaf test = new DendrogramLeaf(5, geneData);
            System.out.println(test);
        }
        catch (Throwable e) {
            System.out.println("Exeption " + e + " caught");
            throw e; // Force improper termination, so error is obvious
        }
    }
}