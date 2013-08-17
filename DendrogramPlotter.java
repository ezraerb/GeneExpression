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

/* This class plots the dendrogram. Its plotted layer by layer from left to
   right starting with the root. The space between leaves in the plot must
   match up with the size of the gene labels in the final layout so the leaves
   and their genes line up correctly */
/* NOTE: Much of this duplicates existing libraries. I wrote it to learn more
   about Java graphics */
public final class DendrogramPlotter extends JPanel
{
    // If the fields in this object change, increment this number by 1
    private static final long serialVersionUID = 2L;

    private DendrogramInterface _dendrogram;
    private int _leafWidth;
        
    /* Constructor. Takes the width between leaves in pixels and the dendrogram
       to plot 
       WARNING: This object takes a reference to the dendrogram. Modifying it
       between constuction and drawing the matrix will lead to invalid
       results! */
    public DendrogramPlotter(DendrogramInterface dendrogram,
                             int leafWidth)
    {
        // Parameter check: both are valid
        if ((dendrogram == null) || (leafWidth <= 0))
            throw new InvalidParameterException();

        _dendrogram = dendrogram;

        _leafWidth = leafWidth;

        /* Set the size. Width is based on number of levels, height on number
           of genes (which is also the number of leaves) */
        Dimension graphSize = new Dimension(dendrogram.getLevel() * _leafWidth,
                                            dendrogram.getGeneCount() * _leafWidth);
        setPreferredSize(graphSize);
        setMinimumSize(graphSize);
    }

    // Draws a dendrogram graphically with the top node at the given position
    private void drawDendrogram(DendrogramInterface dendrogram, Point2D pos,
                                Graphics2D g2)
    {
        if (dendrogram.getLevel() == 0)
            // Leaf node, nothing to do!
            return;

        /* Need to draw lines from the current node to each child. The child
           endpoint X position is one cell size to the right of the current
           node. Finding the Y position is trickier. The final leaves need to
           line up with the gene names, so the position of the node vertically
           depends on how many genes it has relative to the current node.
              The current node contains A genes, of which B came from the first
           child and A - B from the second. The position of the child endpoint
           should be in the exact middle of the genes within the node (remember
           that each leaf lines up with the label for the corresponding gene).
             parent pos = A/2 within the space its part of the correlogram
                         covers
             first child should be B/2 from the top of that space =
                parent pos - A/2 + B/2 = parent - (A - B)/2
             second child should be (A - B)/2 from the bottom of that space =
                parent pos + A/2 - (A - B)/2 = parent pos + A/2 - A/2 + B/2 =
                parentpos + B/2
           Note the surprising appearence of each child node's size in the
           OTHER child's position calculation */
        double firstYadjust = ((double)dendrogram.getSecondChild().getGeneCount() / 2.0) * _leafWidth;
        Point2D firstChildPos = new Point2D.Double(pos.getX() + _leafWidth,
                                                   pos.getY() - firstYadjust);
        double secondYadjust = ((double)dendrogram.getFirstChild().getGeneCount() / 2.0) * _leafWidth;
        Point2D secondChildPos = new Point2D.Double(pos.getX() + _leafWidth,
                                                    pos.getY() + secondYadjust);
        g2.draw(new Line2D.Double(pos, firstChildPos));
        g2.draw(new Line2D.Double(pos, secondChildPos));

        /* Finally, if a child is more than one level below the current node,
           need to draw a horizontal line to cover the skipped levels */
        if (dendrogram.getFirstChild().getLevel() + 1 < dendrogram.getLevel()) {
            Point2D realChildPos = new Point2D.Double(firstChildPos.getX() +
                                                      ((double)(dendrogram.getLevel() - dendrogram.getFirstChild().getLevel() - 1) * _leafWidth),
                                                      firstChildPos.getY());

            g2.draw(new Line2D.Double(firstChildPos, realChildPos));
            // Move the location for drawing the next part of the tree
            firstChildPos = realChildPos;
        } // First child not one level below parent
        else if (dendrogram.getSecondChild().getLevel() + 1 < dendrogram.getLevel()) {
            Point2D realChildPos = new Point2D.Double(secondChildPos.getX() +
                                                      ((double)(dendrogram.getLevel() - dendrogram.getSecondChild().getLevel() - 1) * _leafWidth),
                                                      secondChildPos.getY());

            g2.draw(new Line2D.Double(secondChildPos, realChildPos));
            // Move the location for drawing the next part of the tree
            secondChildPos = realChildPos;
        } // Second child not one level below parent

        // Draw the levels below the children
        drawDendrogram(dendrogram.getFirstChild(), firstChildPos, g2);
        drawDendrogram(dendrogram.getSecondChild(), secondChildPos, g2);
    }

    public void paintComponent(Graphics g)
    {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D)g;

        /* Dendrograms are trees, which are best processed recursively. Set up
           the overall parameters, and then start on the dendrogram root. */
        Stroke saveStroke = g2.getStroke();
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                            RenderingHints.VALUE_ANTIALIAS_ON);
        g2.setStroke(new BasicStroke(3));
        
        /* The tree starts on the left edge, halfway down the panel, which is
           set from the number of genes in the tree */
        Point2D root = new Point2D.Double(0.0,
                                          _leafWidth * _dendrogram.getGeneCount() / 2.0);
        drawDendrogram(_dendrogram, root, g2);
        g2.setStroke(saveStroke); // Restore original so rest of graph looks right
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
            DendrogramPlotter test = new DendrogramPlotter(cluster.getDendrogram(),
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
