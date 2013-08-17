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

/* These interface defines a dendrogram, a tree that groups genes by similarity
   of expression pattern. Its implemented using the classic whole-part design
   patern, so all nodes act identically. This allows the dendrogram to be cut
   on any level to get the grouping of genes as of that level. The nodes also
   have a sort order based on the expression data of the genes within the node;
   a full traversal of the tree will produce a list of genes sorted by
   expression similarity */
public interface DendrogramInterface
{
    /* Number of genes in this node or its children. A value of one indicates
       a leaf */
    public int getGeneCount();

    // Names of genes in the node
    public Set<String> getGenes();
        
    /* Total of the weigts of all leaves under this node, or the weight if
       its a leaf */
    public double getTotalWeight();

    /* Averegate weight of the node, defined as the average leaf weight, used
       for sorting */
    public double getAverageWeight();

    /* Level of this node, defined as the number of nodes to reach the
       furthest leaf */
    public int getLevel();
    
    // Returns first child, if any
    public DendrogramInterface getFirstChild();
    
    // Returns second child, if any
    public DendrogramInterface getSecondChild();

    // Output the portion of the tree of this node and below
    public void toString(StringBuffer output, int maxLevel, int parentLevel);
    
    // Output the tree starting with this node
    public String toString();
}
