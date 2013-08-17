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

public final class GeneExpression extends JFrame
{
    // Constructor. Needs the clustered data
    public GeneExpression(Cluster data, Dimension graphSize)
    {
        if (data == null)
            throw new InvalidParameterException("No gene data passed to graph");

        // Extract the clustered corrolation matrix, needed to build the graph
        GeneCorrolations geneData = data.getSortedCorrolations();
        DendrogramInterface dendrogram = data.getDendrogram(); 

        // Sanity check
        if (geneData.getGeneCount() == 0)
            throw new InvalidParameterException("No gene data passed to graph");
                    
        /* The final graph has three components, the denndrogram tree on the
           left, the gene labels in the middle, and the correlogram on the
           right. To ensure everything lines up, need the distance between
           gene names, which depends on how they will render. Find the longest
           gene name, construct a label with it, and get the size of the text
           given the font. Thankfully, the label does not need to actually be
           rendered to get the text size.
           WARNING: This assumes that letter sizes in fonts are approximately
           equal */
        String longestGene = geneData.getGeneName(0); // Guarenteed to exist
        int index;
        for (index = 1; index < geneData.getGeneCount(); index++) {
            String testGene = geneData.getGeneName(index);
            if (longestGene.length() < testGene.length())
                longestGene = testGene;
        } // Loop through gene names

        JLabel testLabel = new JLabel(longestGene);
        FontMetrics testMetrics = testLabel.getFontMetrics(testLabel.getFont());

        int labelWidth = testMetrics.stringWidth(longestGene);
        int labelHeight = testMetrics.getHeight();
        
        /* Increase both by 10% to ensure a margin of error. The funky math is
           due to making integers round properly */
        labelWidth = ((labelWidth * 11) + 5) / 10;
        labelHeight = ((labelHeight * 11) + 5) / 10;

        /* Construct the graph using grid-bag layout. The components are sized
           based on the hight of the gene name labels, ensuring everything
           lines up correctly */
        JPanel finalGraph = new JPanel(new GridBagLayout());
        GridBagConstraints graphLayout = new GridBagConstraints();
        graphLayout.gridy = 0;
        graphLayout.gridx = 0;
        graphLayout.gridheight = geneData.getGeneCount();
        finalGraph.add(new DendrogramPlotter(dendrogram, labelHeight),
                       graphLayout);
        
        graphLayout.gridx = 2;
        finalGraph.add(new CorrolationPlotter(geneData, dendrogram,
                                              labelHeight),
                       graphLayout);
        

        // Iterate through the gene names and generate the labels
        Dimension labelSize = new Dimension(labelWidth, labelHeight);
        graphLayout.gridx = 1;
        graphLayout.gridheight = 1;
        for (index = 0; index < geneData.getGeneCount(); index++) {
            graphLayout.gridy = index;
            JLabel newLabel = new JLabel(geneData.getGeneName(index),
                                         SwingConstants.CENTER);
            newLabel.setSize(labelSize);
            newLabel.setPreferredSize(labelSize);
            finalGraph.add(newLabel, graphLayout);
        } // Loop through gene names

        /* Most gene corrolations graphs are so large they can not possibly
           fit on the screen. Use a scroll pane to ensure everything displays */
        getContentPane().add(new JScrollPane(finalGraph), BorderLayout.CENTER);
        setTitle("Gene Corrolations");
        setSize(graphSize);
        pack();
        setVisible(true);
    }

    /* Find the maximum size available for the plot window. This is the screen
       size minus taskbars and the like */
    private static Dimension getScreenSize()
    {
        GraphicsConfiguration screenData = GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice().getDefaultConfiguration();
        Rectangle screenSize = new Rectangle(screenData.getBounds());
        Insets insets = Toolkit.getDefaultToolkit().getScreenInsets(screenData);
        screenSize.height -= (insets.top + insets.bottom);
        screenSize.width -= (insets.left + insets.right);

        /* Rectangle uses double values while the frame sizing routines want
           integer values. Hence this clunky implementation */
        return new Dimension((int)screenSize.width, (int)screenSize.height);
    }

    /* Main program. Reads a data file from the command line, clusters the
       data it contains, and displays the results */
    public static void main(String[] args) throws Throwable
    {
        try {
            if (args.length != 1)
                throw new IllegalArgumentException("Wrong number of arguments, need file name only");
            DataLoader data = new DataLoader();
            data.loadData("data", args[0], 5.0);
            Cluster cluster = new Cluster(new GeneCorrolations(data));
            GeneExpression output = new GeneExpression(cluster,
                                                       getScreenSize());
            new JFrameThreadWrapper(output, 1);
        }
        catch (Throwable e) {
            System.out.println("Exeption " + e + " caught");
            throw e; // Force improper termination, so error is obvious
        }
    }
}