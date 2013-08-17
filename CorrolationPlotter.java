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
import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.util.*;
import java.security.InvalidParameterException;

/* This class plots the corrolations maxtrix. The shading of a gene pair
   represents their level of expression similiarity. Color blocks indicate
   areas of greatest similiarity based on the dendrogram. The size of cells
   much match up with the gene labels so the rows are labeled properly */
/* NOTE: Much of this duplicates existing libraries. I wrote it to learn more
   about Java graphics */
public final class CorrolationPlotter extends JPanel
{
    // If the fields in this object change, increment this number by 1
    private static final long serialVersionUID = 2L;

    private int _cellSize;
    private GeneCorrolations _corrolations;
    
    /* Enumerated type to represent wanted colors. This is a deleberate subset
       of all available restricted to those which look nice when shaded */
    private enum DisplayColors
    {
        BLUE(Color.blue),
        CYAN(Color.cyan),
        GREEN(Color.green),
        MAGENTA(Color.magenta),
        ORANGE(Color.orange),
        PINK(Color.pink),
        RED(Color.red),
        YELLOW(Color.yellow);

        Color _color; // Color this enum value represents

        // Constructor
        DisplayColors(Color color)
        {
            _color = color;
        }

        public Color getColor()
        {
            return _color;
        }
    };

    // Set of all of the above. Needed to iterate it later
    /* NOTE: Java Enums have no built in method to iterate through the values,
       other than using a for loop. That mechanism won't meet this code's
       needs, because the iteration happens based on gene groups. Using a set
       and an iterator is a standard design pattern (some call it an
       anti-pattern) for handling the situation */
    private Set<DisplayColors> _colorList;
        
    /* List of indexes of genes in the corrolation matrix considered to be
       within a single output group (found from the dendrogram) */
    private int[] _geneGroupLimits;

    /* Constructor. Takes the width of each cell in pixels, the sorted and
       normalized corrolation matrix, and the corresponding dendrogram
       WARNING: This object takes a reference to the corrolation matrix.
       Modifying it between constuction and drawing the matrix will lead to
       invalid results! */
    public CorrolationPlotter(GeneCorrolations corrolations,
                              DendrogramInterface dendrogram, int cellSize)
    {
        // Parameter check: all are valid
        if ((corrolations == null) || (dendrogram == null) || (cellSize <= 0))
            throw new InvalidParameterException();
        /* Consistency check: number of genes in both corrolation matrix and
           dendrogram match */
        if (corrolations.getGeneCount() != dendrogram.getGeneCount())
            throw new InvalidParameterException();

        _corrolations = corrolations;

        _colorList = EnumSet.allOf(DisplayColors.class);

        /* Group genes by similarity from the dendrogram. Set the limit to the
           lesser of the number of group colors and half the number of genes */
        int groupLimit = dendrogram.getGeneCount() / 2;
        if (groupLimit > _colorList.size())
            groupLimit = _colorList.size();
        _geneGroupLimits = groupGenes(dendrogram, groupLimit);
        
        _cellSize = cellSize;

        // Set the size. Width and height is based on the number of genes
        Dimension graphSize = new Dimension(corrolations.getGeneCount() * _cellSize,
                                            corrolations.getGeneCount() * _cellSize);
        setPreferredSize(graphSize);
        setMinimumSize(graphSize);
    }

    /* Collect genes into groups for output. The most similiar genes, based
       on the correlogram, should be colored the same in the output. */
    private int[] groupGenes(DendrogramInterface dendrogram,
                             int groupCount)
    {
        /* Finding proximate gene groups by a dendrogram is based on a few
           conceps:
           1. Every node of the dendrogram (including the root) represents a
           grouping of genes as of a level of proximity
           2. The set of genes represented by a given node is the combination
           of the sets of genes of its children.
           3. The higher the node level, the lower the proximity amount of the
           nodes within it
           They mean that finding the gene groups is the same as finding a set
           of nodes from the dendrogram; such that no node is above any other
           node in the tree, the levels of the nodes are as low as possible,
           and their genes collectively cover all in the dendrogram.

           The way this method finds it is by starting with a list of only one
           node, the root. It contains all genes. The root is then replaced in
           the list by its children. This effectively splits the gene group
           into two. One of the children is selected and replaced with its
           children, creating three groups. The process repeats until the
           wanted number of nodes appears in the set.

           Choosing which node to split is critical. The obvious choice is the
           node with the highest level. After a few splits, multiple nodes will
           likely have the same level. Ties are broken by taking the node with
           lower level children, then the node with more genes in it.

           Data structure discussion: The group sort order in the final results
           is the same as the order of the nodes in a tree transversal. Sorting
           them at the end will be tough, so this algorithm produces them in
           order. The children of a replaced node are inserted in the node
           list in the same location as their parent, with the first child
           first. A linked list allows constant insert at the cost of a linear
           search. The number of groups will usually be low enough the compute
           is tolerable.

           The genes will be listed in the corrolelogram in the same order as
           the leaves of the dendrogram. This means the final gene groups can
           be converted to index values of the last gene in the group to group
           correlelogram entries for graphing */
        LinkedList<DendrogramInterface> nodeList = new LinkedList<DendrogramInterface>();
        nodeList.add(dendrogram);

        /* NOTE: Not worth checking for a group count larger than then the
           number of genes, it won't happen often in practice */
        int nextNode = 0;
        while ((nodeList.size() < groupCount) && (nextNode != -1)) {
            // Find next node to split to create another group
            nextNode = -1;
            // Find next node to replace. This requires a linear scan
            int childLevel = 0;
            int geneCount = 0;
            int testIndex;
            for (testIndex = 0; testIndex < nodeList.size(); testIndex++) {
                DendrogramInterface testNode = nodeList.get(testIndex);
                if (testNode.getGeneCount() > 1) { // Branch node
                    int testChildLevel = testNode.getFirstChild().getLevel();
                    if (testChildLevel < testNode.getSecondChild().getLevel())
                        testChildLevel = testNode.getSecondChild().getLevel();
                    if ((nextNode == -1) || // No replace node found yet
                        (testNode.getLevel() > nodeList.get(nextNode).getLevel()) ||
                        ((testNode.getLevel() == nodeList.get(nextNode).getLevel()) &&
                         ((testChildLevel < childLevel) ||
                          ((testChildLevel == childLevel) &&
                           (testNode.getGeneCount() > geneCount))))) {
                        nextNode = testIndex;
                        childLevel = testChildLevel;
                        geneCount = testNode.getGeneCount();
                    } // This node is the best found so far
                } // For each branch node
            } // For each node in the current list

            /* If a branch to process was found, replace it with its children.
               Due to how linked lists work, overwrite the current node with
               its second child, and then insert the first */
            if (nextNode != -1) {
                DendrogramInterface tempNode = nodeList.get(nextNode);
                nodeList.set(nextNode, tempNode.getSecondChild());
                nodeList.add(nextNode, tempNode.getFirstChild());
            } // Node to update found
        } // While nodes to replace and reasons to do so

        /* Convert to index values. The dendrogram groups are in the same
           order as the correlogram rows, so the size of each group gives
           the highest index in the group */
        int[] results = new int[nodeList.size()];
        int currIndex = 0;
        for (nextNode = 0; nextNode < nodeList.size(); nextNode++) {
            currIndex += nodeList.get(nextNode).getGeneCount();
            results[nextNode] = currIndex;
        }
        return results;
    }

    /* Given a color and a darkness value, returns a new color darkended by
       the wanted amount */
    private Color colorBrightness(Color original, float brightness)
    {
        /* Surprisingly, Java has no method to apply a brightness to a default
           RGB color. This method provides it. */

        // Ensure consistency regardless of original color model
        float [] components = new float[3];
        components[0] = (float)original.getRed() / 255;
        components[1] = (float)original.getGreen() / 255;
        components[2] = (float)original.getBlue() / 255;
        int index;
        for (index = 0; index < components.length; index++)
            components[index] *= brightness;
        return new Color(components[0], components[1], components[2]);
    }
    
    // Draws one cell of the graph
    private void drawCell(int row, int column, Graphics2D g2)
    {
        /* Yes, row and column look swapped. The rectable wants the X and Y
           cordinates of the upper right corner, which depend on the column
           and row, respectively */
        g2.fill(new Rectangle2D.Double((double)column * _cellSize,
                                       (double)row * _cellSize,
                                       (double)_cellSize,
                                       (double)_cellSize));
    }

    public void paintComponent(Graphics g)
    {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D)g;

        /* Get an iterator to the first entry of the set of colors to use.
           When the row index points to a gene outside the current set of
           genes to group from the dendrogram, increment the color
           pointer and move to the next group */
        Iterator<DisplayColors> colorIndex = _colorList.iterator();
        DisplayColors groupColor = colorIndex.next();
        int currGroup = 0;
        int firstGroupRow = 0;
        
        /* Iterate through the corrolations matrix graph and output each cell
           with the brightness based on the cell value. A LOWER value indicates
           a higher corrolation, so subtract from 1 before applying it. The
           graph is symetric, so output rows and columns simultaneously */
        int rowIndex, colIndex;
        for (rowIndex = 0; rowIndex < _corrolations.getGeneCount();
             rowIndex++) {
            if (rowIndex >= _geneGroupLimits[currGroup]) {
                // Move to next group
                groupColor = colorIndex.next();
                // End of old group is start of next
                firstGroupRow = _geneGroupLimits[currGroup]; 
                currGroup++;
            }
            /* Given the order squares are drawn, the group squares, closest to
               the diagonal, are always drawn last. Start the color white (for
               squares outside the group) and then switch it at the appropriate
               column */
            Color cellColor = Color.white; 
            for (colIndex = 0; colIndex < rowIndex; colIndex++) {
                if (colIndex >= firstGroupRow)
                    cellColor = groupColor.getColor();
                g2.setPaint(colorBrightness(cellColor, 
                                            ((float)(1.0 - _corrolations.getCorrolation(rowIndex, colIndex)))));
                drawCell(rowIndex, colIndex, g2);
                drawCell(colIndex, rowIndex, g2);
            } // Column loop
            // A gene is always fully corrolated with itself, and inside its group
            g2.setPaint(groupColor.getColor());
            drawCell(rowIndex, rowIndex, g2);
        } // loop through rows
    }

    // Object test program. Reads a hard coded data file and clusters it
    public static void main(String[] args) throws Throwable
    {
        try {
            DataLoader data = new DataLoader();
            data.loadData("data", "test.txt", 5.0);
            GeneCorrolations geneData = new GeneCorrolations(data);
            Cluster cluster = new Cluster(geneData);

            // Create the corrolation plot with the results
            CorrolationPlotter test = new CorrolationPlotter(cluster.getSortedCorrolations(),
                                                             cluster.getDendrogram(),
                                                             18);
            GraphTestFramework test2 = new GraphTestFramework(test, 500, 500);
            test2.setVisible(true);
        }
        catch (Throwable e) {
            System.out.println("Exeption " + e + " caught");
            throw e; // Force improper termination, so error is obvious
        }
    }
}
    