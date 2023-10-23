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
package application;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class InputParametersParser {

	static public List<String> input = new ArrayList<String>();
	static public List<String> bootstrap = new ArrayList<String>();
	static public List<String> target = new ArrayList<String>();
	static public List<String> userDefinedPrecursor = new ArrayList<String>();
	static public List<String> forbiddenPrecursors = new ArrayList<String>();
	
	static public void parse(String inputFile) {
		input.clear();
		bootstrap.clear();
		target.clear();
		userDefinedPrecursor.clear();
		
		boolean inputs = false;
		boolean bootstraps = false;
		boolean targets = false;
		boolean udPrecursors = false;
		boolean forbiddenCompounds = false;
		
		String line = "";
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(inputFile));
			line = br.readLine().trim();
			
			while(line != null)
			{
				line = line.trim();
				if( line.equalsIgnoreCase("<input-compounds>") )
				{
					inputs = true;
					bootstraps = false; targets = false; udPrecursors = false; forbiddenCompounds = false;
				}
				else if( line.equalsIgnoreCase("<bootstrap-compounds>") )
				{
					bootstraps = true;
					inputs = false; targets = false; udPrecursors = false; forbiddenCompounds = false;
				}
				else if( line.equalsIgnoreCase("<precursor-compounds>") )
				{
					udPrecursors = true;
					inputs = false; bootstraps = false; targets = false; forbiddenCompounds = false;
				}
				else if( line.equalsIgnoreCase("<target-compounds>") )
				{
					targets = true;
					inputs = false; bootstraps = false; udPrecursors = false; forbiddenCompounds = false;
				}
				else if( line.equalsIgnoreCase("<forbidden-compounds>")) 
				{
					forbiddenCompounds = true;
					inputs = false; bootstraps = false; udPrecursors = false; targets = false;
				}
				else if( line.startsWith("<species") ) {
					String startingWithId = line.substring(line.indexOf("=")+2);
					String id = startingWithId.substring(0, startingWithId.indexOf('"'));
					if( inputs )
						input.add(id);
					else if (bootstraps )
						bootstrap.add(id);
					else if( udPrecursors )
						userDefinedPrecursor.add(id);
					else if( targets )
						target.add(id);
					else if( forbiddenCompounds )
						forbiddenPrecursors.add(id);
				}
				line = br.readLine();
			}
			br.close();
		} catch (FileNotFoundException e) {
			System.out.println("Input file not found!");
			e.printStackTrace();
			System.exit(-1);
		} catch (IOException e) {
			System.out.println("Error reading input!");
			System.out.println("Last line: "+line);
			e.printStackTrace();
			System.exit(-1);
		}
		
	}

	public static List<String> getInputForFP() {
		return input;
	}

	public static void setInputForFP(List<String> inputForFP) {
		InputParametersParser.input = inputForFP;
	}

	public static List<String> getBootstrapForFP() {
		return bootstrap;
	}

	public static void setBootstrapForFP(List<String> bootstrapForFP) {
		InputParametersParser.bootstrap = bootstrapForFP;
	}

	public static List<String> getTarget() {
		return target;
	}

	public static void setTarget(List<String> target) {
		InputParametersParser.target = target;
	}

	public static List<String> getUserDefinedPrecursor() {
		return userDefinedPrecursor;
	}

	public static void setUserDefinedPrecursor(List<String> userDefinedPrecursor) {
		InputParametersParser.userDefinedPrecursor = userDefinedPrecursor;
	}

	public static List<String> getForbiddenPrecursors() {
		return forbiddenPrecursors;
	}

	public static void setForbiddenPrecursors(List<String> forbiddenPrecursors) {
		InputParametersParser.forbiddenPrecursors = forbiddenPrecursors;
	}
	
}
