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
import java.io.*;
import java.util.*;
import java.lang.*;

/* This class loads gene expression data from a data file and normalizes it.
   Minor errors and inconsisties in the data are flagged and ignored. Severe
   errors will cause termination */
public final class DataLoader {

    // File to load plays. Inside class to ensure always released
    private BufferedReader _file;

    /* Data from the file. The final number of entries is unknown, hence
       the need for a dynamic data structure. */
    private ArrayList<GeneData> _geneData;

    public DataLoader()
    {
        // Construct object empty
        _geneData = new ArrayList<GeneData>();
    }

    /* Processes a line of data from the file and updates results data.
       Returns status of processing */
    private boolean processFileLine(String buffer, String fullName,
                                    int dataFileLine, double tolerance,
                                    boolean noErrorsSoFar)
    {
        /* SUBTLE NOTE: Want to process the file line even if other errors
           have already been found, in order to log as many as possible. */
        boolean success = true;
        String [] newData = buffer.split(" ");
        if (newData.length <= 1) {
            System.out.println("File " + fullName + " corrupt; line " +
                               dataFileLine + ": " + buffer +
                               ": observations missing");
            success = false;
        }
        else {
            // Convert the observations as strings to numeric values
            double [] newResult = new double[newData.length - 1];
            int index;
            for (index = 0; index < newResult.length; index++) {
                try {
                    newResult[index] = Double.parseDouble(newData[index + 1]);
                }
                catch (Exception e) {
                    // Conversion error. Log and flag
                    System.out.println("File " + fullName + " corrupt; line " +
                                       
                                       dataFileLine + ": " + buffer +
                                       " observation " + index + " invalid");
                    success = false;
                    newResult[index] = 0.0; // Ensure consistency
                } // Conversion error
            } // For loop

            /* If no errors in this line, confirm the number of observations is
               consistent with the number for previous genes */
            if (success && (!_geneData.isEmpty())) {
                if (_geneData.get(0).getNumObservations() != newResult.length) {
                    System.out.println("File " + fullName + " corrupt; line " +
                                       
                                       dataFileLine + ": " + buffer +
                                       " expected observations: " +
                                       _geneData.get(0).getNumObservations() +
                                       " actual: " + newResult.length);
                    success = false;
                }
            } // Previous observations

            /* If no errors in this line, or the previous file lines, normalize
               the results. If the resulting variation is below the passed
               tolerance, drop the gene because it's expression did not vary
               sufficiently */
            /* SUBTLE NOTE: If the previous file lines had errors, don't
               bother doing the insert because the file results will ultimately
               be discarded */
            if (success && noErrorsSoFar) {
                double variance = normalizeData(newResult);
                if (variance >= tolerance)
                    _geneData.add(new GeneData(newData[0], newResult));
            }
        } // Line read from file contains observations

        // Return the status for the overall file
        return (noErrorsSoFar && success);
    }

    /* Normalizes a set of observations and returns the variance of the original
       observations. Its calculated as part of the normalization process, so
       just return it instead of trying to calculate it twice */
    private double normalizeData(double [] data)
    {
        /* Gene expression calculations care about the relative expression
           levels of different genes. Using the raw values distorts the
           results, so the data is almost always normalized first. The usual
           algorithm is to normalize to a mean of zero and a variance of one,
           which mimics a well used correlation formula called Pearson's
           Coefficient. */
        if (data.length == 0)
            return 0.0; // Avoid divide by zero error

        /* First, normalize the distribution to mean zero. Need the current mean
           to do so. Note that this normalization does not affect the standard
           deviation value */
        double temp = 0.0; // Higher precision to avoid overflow
        int index;
        for (index = 0; index < data.length; index++)
            temp += data[index];
        double mean = temp / data.length;
        for (index = 0; index < data.length; index++)
            data[index] -= mean;

        /* Find the standard deviation. Note that the adjustment above did not
           change the value. The calculation is faster on a zero-mean
           distribution, so done here */
        temp = 0.0;
        for (index = 0; index < data.length; index++)
            temp += (data[index] * data[index]);
        double deviation = Math.sqrt(temp / data.length);

        /* To normalize a distribution centered on zero, divide by the
           standard deviation */
        if (deviation != 0.0)
            for (index = 0; index < data.length; index++)
                data[index] /= deviation;
        return deviation;
    }

    /* Loads data from the passed file into the object. Duplicate gene data
       is not checked for. If the number of observations per gene is not
       consistent, or multiple files are loaded and the number of observations
       is not consistent across files, an exception is issued. Genes whose
       expression levels have insignificant variation are automatically
       discarded. Observations are normalized */
    public void loadData(String directory, String fileName, double tolerance)
        throws Exception
    {
        // Record current size of vectors, so location of new data is known
        int currDataSize = _geneData.size();
        try {
            // Open the file. Not finding it is an error
            String fullName;
            if (directory != null) // Directory passed
                /* TRICKY NOTE: Notice the double backslash below. Java uses
                   '\' as an escape character. The first is the esacpe
                   character needed to insert a litteral '\' in the string! */
                fullName = directory + '\\' + fileName;
            else
                fullName = fileName;
            _file = new BufferedReader(new FileReader(fullName));
            if (_file == null) {
                // Failed to open. This is a big error 
                System.out.println("Error loading play data, file " + fullName + " missing");
                throw new IOException("Error loading play data, file " + fullName + " missing");
            } // File not opened
            else {
                // Read all lines even for errors so they are all logged
                boolean success = true;

                /* Single lines of data may break over multiple file lines,
                   so need to recombine them to process them */
                String fileData = null;

                // First line is a header. Read it to burn it
                _file.readLine();
                String buffer = _file.readLine();

                // Track lines containing data to report locations in errors
                int fileLine = 1;
                int dataFileLine = 0;
                while (buffer != null) {

                    /* If the buffer contains only spaces (or nothing at all)
                       have the end of a data line. If read something,
                       process it */
                    if (buffer.trim().isEmpty()) {
                        if (fileData != null) {
                            success = processFileLine(fileData, fullName,
                                                      dataFileLine, tolerance,
                                                      success);
                            fileData = null;
                        }
                        // Toss the empty line
                    }
                    else {
                        // Data. If no file line found yet, have it now
                        if (fileData == null) {
                            fileData = buffer;
                            dataFileLine = fileLine;
                        }
                        /* If the file line contains things other than numbers
                           and spaces, have a new set of observations. Process
                           the existing set, and then copy the new data over */
                        else if (buffer.matches("[^\\d\\s]+")) {
                            success = processFileLine(fileData, fullName,
                                                      dataFileLine, tolerance,
                                                      success);
                            fileData = buffer;
                            dataFileLine = fileLine;
                        }
                        else
                            /* More of the same set of observations. Append
                               to the previous set */
                            fileData = fileData.concat(" " + buffer);
                    } // Non-empty line from file
                    buffer = _file.readLine();
                    fileLine++;
                } // While lines in the file to process

                /* May get to here without having processed the last line of
                   data. Handle it now */
                if (fileData != null) {
                    success = processFileLine(fileData, fullName, dataFileLine,
                                              tolerance, success);
                    fileData = null;
                }

                // If have a data error, clear the new data from the results
                if ((!success) && (currDataSize < _geneData.size()))
                    _geneData.subList(currDataSize, _geneData.size()).clear();
                closeFile();
            } // File successfully opened
        } // Try..catch block around entire routine
        catch (Exception e) {
            /* Clean up before continuing */
            closeFile();
            // If any entries were inserted into the data, remove them
            if (currDataSize < _geneData.size())
                _geneData.subList(currDataSize, _geneData.size()).clear();
            throw e;
        }
    }

    /* Closes the file open in the object */
    private void closeFile() throws Exception
    {
        if (_file != null)
            _file.close();
        _file = null;
    }

    /** This method ensures the file is always closed before the object dies.
        In general, if the file gets to here, something has gone wrong and
        resources have been held far longer than needed. A warning is issued
        to handle this case */
    public void finailize()
    {
        if (_file != null) {
            System.out.println("WARNING: File not properly closed");
            try {
                _file.close();
            }
            /* Catch and dispose of any IO exception, since the object is going
               away soon anyway. This is normally an anti-pattern, but needed
               in this case */
            catch (IOException e) {}
            _file = null;
        }
    }

    // Dumps the contects of the class to a string
    public String toString()
    {
        if (_geneData.isEmpty())
            return ("no data loaded");
        else {
            StringBuffer output = new StringBuffer();
            Iterator<GeneData> index = _geneData.iterator();
            while (index.hasNext()) {
                output.append(index.next());
                output.append('\n');
            }
            return output.toString();
        } // Data loaded
    }

    // Returns the data within the object. WARNING the actual data, not a copy
    public ArrayList<GeneData> getGeneData()
    {
        return _geneData;
    }

    // Object test program. Reads a hard coded data file, outputs the results
    public static void main(String[] args) throws Throwable
    {
        try {
            DataLoader test = new DataLoader();
            test.loadData("data", "test.txt", 5.0);
            System.out.println(test);
        }
        catch (Throwable e) {
            System.out.println("Exeption " + e + " caught");
            throw e; // Force improper termination, so error is obvious
        }
    }
}