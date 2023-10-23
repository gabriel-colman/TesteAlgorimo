/* 
PITUFO - A software tool to find all minimal precursor sets for a given set of targets in metabolic networks.

Copyright (C) 2011 Ludovic Cottret (l.cottret@gmail.com <mailto:l.cottret@gmail.com>), Paulo Vieira Milreu  (paulovieira@milreu.com.br <mailto:paulovieira@milreu.com.br>)      
This file is part of PITUFO.

PITUFO is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PITUFO is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PITUFO.  If not, see <http://www.gnu.org/licenses/>.
*/
package utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class FileUtils {

	// read a file and store the content in a list (array element == a line)
	public static List<String> readTextFile(String filename){
		List<String> content = new ArrayList<String>();
		BufferedReader br = null;
		 
		try {
 
			String sCurrentLine;
 
			br = new BufferedReader(new FileReader(filename));
 
			while ((sCurrentLine = br.readLine()) != null) {
				content.add(sCurrentLine);
			}
 
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
		return content;
	}
	
	
	  public static void emptyDirectory(String directoryName) 
	  {
	    File dir = new File(directoryName);
	    if (!dir.exists()) {
	      System.out.println(directoryName + " does not exist");
	      return;
	    }
	
	    String[] info = dir.list();
	    for (int i = 0; i < info.length; i++) 
	    {
	    	File n = new File(directoryName + File.separator + info[i]);
	    	if (!n.isFile()) // skip ., .., other directories too
	    		continue;
	    	System.out.println("removing " + n.getPath());
	    	if (!n.delete())
	    		System.err.println("Couldn't remove " + n.getPath());
	    }
	  }
	  
	  public static void eraseFilesByName(String directoryName, String pattern) 
	  {
	    File dir = new File(directoryName);
	    if (!dir.exists()) {
	      System.out.println(directoryName + " does not exist");
	      return;
	    }
	    String[] info = dir.list();
	    for (int i = 0; i < info.length; i++) {
	      File n = new File(directoryName + File.separator + info[i]);
	      if (!n.isFile()) { // skip ., .., other directories, etc.
	        continue;
	      }
	      if (info[i].indexOf(pattern) == -1) { // name doesn't match
	        continue;
	      }
	      System.out.println("removing " + n.getPath());
	      if (!n.delete())
	        System.err.println("Couldn't remove " + n.getPath());
	    }
	 }	  
}
