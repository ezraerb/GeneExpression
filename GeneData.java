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

// This class holds normalized gene expression observation data
public final class GeneData
{
    private String _geneName;
    private double [] _geneData;

    // Constructor
    public GeneData(String geneName, double [] geneData)
    {
        _geneName = geneName;
        _geneData = geneData;
    }
    
    // Getters
    public String getGeneName()
    {
        return _geneName;
    }
    public double [] getGeneData()
    {
        return _geneData;
    }
    
    // Returns number of observations
    public int getNumObservations()
    {
        return _geneData.length;
    }
    
    // Dumps data to a string
    public String toString()
    {
        return _geneName + ":" + Arrays.toString(_geneData);
    }
}

