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

/* This class defines a branch of a dendrogram, a tree that groups genes by
   similarity of expression pattern. Its implemented using the classic
   whole-part design patern, so all nodes act identically. This allows the
   dendrogram to be cut on any level to get the grouping of genes as of that
   level. The nodes also have a sort order based on the expression data of the
   genes within the node */
public final class DendrogramBranch implements DendrogramInterface
{
    private DendrogramInterface _firstChild;
    private DendrogramInterface _secondChild;

    /* Data for this node is calculated by combining the children. Since that
       calculation is recursive, recalculating it for every method call will
       get expensive fast. Precalculate it and cache it here instead. Children
       are only set in the constructor, so don't have to worry about cache
       updates */
    private Set<String> _genes;

    // Node weight, defined as the total of the weights of child nodes
    private double _weight;

    // Level. One higher than the highest level of the child nodes
    private int _level;

    // Constructor for a dendrogram branch. It requires two children 
    public DendrogramBranch(DendrogramInterface firstChild,
                            DendrogramInterface secondChild)
    {
        if ((firstChild == null) || (secondChild == null))
            throw new InvalidParameterException();

        /* Construct the set of all genes in this branch node, which is the
           combination of all genes in its children. If the size of this set is
           less than the sum of the two children sets, it indicates a gene
           appears in both, which in turn signals a malformed dendrogram
           tree. It can also be caused by bad data reporting a gene twice */
        _genes = new HashSet<String>(firstChild.getGenes());

        _genes.addAll(secondChild.getGenes());
        if (_genes.size() !=
            (firstChild.getGeneCount() + secondChild.getGeneCount())) {

            /* Find the duplicate gene name by taking the intersection of the
               original sets */
            _genes = new HashSet<String>(firstChild.getGenes());
            _genes.retainAll(secondChild.getGenes());
            System.out.println("ERROR: Dendrogram child gene records duplicate " + _genes);
            throw new InvalidParameterException();
        }
            
        _weight = firstChild.getTotalWeight() + secondChild.getTotalWeight();
        
        _level = firstChild.getLevel();
        if (secondChild.getLevel() > _level)
            _level = secondChild.getLevel();
        // Branch is one level higher than children
        _level++;

        /* To get the sort order correct, the first child must be the one with
           the lower average weight */
        if (firstChild.getAverageWeight() <= secondChild.getAverageWeight()) {
            _firstChild = firstChild;
            _secondChild = secondChild;
        }
        
        else {
            _firstChild = secondChild;
            _secondChild = firstChild;
        }
    }

    /* Number of genes in this node or its children. A value of one indicates
       a leaf */
    public int getGeneCount()
    {
        return _genes.size();
    }

    // Names of genes in the node
    public Set<String> getGenes()
    {
        return _genes;
    }
        
    /* Total of the weigts of all leaves under this node, or the weight if
       its a leaf */
    public double getTotalWeight()
    {
        return _weight;
    }

    /* Averegate weight of the node, defined as the average leaf weight, used
       for sorting */
    public double getAverageWeight()
    {
        return _weight / getGeneCount();
    }

    /* Level of this node, defined as the number of nodes to reach the
       furthest leaf */
    public int getLevel()
    {
        return _level;
    }
    
    // Returns first child
    public DendrogramInterface getFirstChild()
    {
        return _firstChild;
    }
    
    // Returns second child
    public DendrogramInterface getSecondChild()
    {
        return _secondChild;
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
        output.append("BRANCH level: " + _level + " Total weight: " + _weight +
                      '\n');
        // Output both children
        _firstChild.toString(output, maxLevel, getLevel());
        _secondChild.toString(output, maxLevel, getLevel());
    }
    
    // Output the tree starting with this node
    public String toString()
    {
        StringBuffer output = new StringBuffer();
        toString(output, getLevel(), 0);
        return output.toString();
    }

    // Object test program. Reads a hard coded data file and manipulates it
    public static void main(String[] args) throws Throwable
    {
        try {
            DataLoader data = new DataLoader();
            data.loadData("data", "test.txt", 5.0);
            GeneCorrolations geneData = new GeneCorrolations(data);
            DendrogramLeaf leaf1 = new DendrogramLeaf(5, geneData);
            DendrogramLeaf leaf2 = new DendrogramLeaf(7, geneData);
            DendrogramLeaf leaf3 = new DendrogramLeaf(0, geneData);
            DendrogramBranch test1 = new DendrogramBranch(leaf1, leaf3);
            DendrogramBranch test2 = new DendrogramBranch(test1, leaf2);
            System.out.println(test2);
        }
        catch (Throwable e) {
            System.out.println("Exeption " + e + " caught");
            throw e; // Force improper termination, so error is obvious
        }
    }
}
