package application;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import metabolicNetwork.Compound;

public class PrecursorSetSolutionParser {

	//static public List<PrecursorSet> precursorSetSolutions = new ArrayList<PrecursorSet>();
	
	static public List<PrecursorSet> parse(File directory, String target) {
		
		List<PrecursorSet> precursorSetSolutions = new ArrayList<PrecursorSet>();
		String inputFile = directory.getPath() + "/solutions_" + target + "_PS.xml";
		try {
			String line = "";
			BufferedReader br = new BufferedReader(new FileReader(inputFile));
			line = br.readLine().trim();
			
			PrecursorSet ps = new PrecursorSet(0, 0, 0);
			while(line != null){
				line = line.trim();
				if(line.startsWith("<source id")){ // create source compound
					String id = line.replaceAll("(.*id=\")(.*?)(\"\\s.*)", "$2");
					String name = line.replaceAll("(.*name=\")(.*?)(\"\\s.*)", "$2");
					String compartment = line.replaceAll("(.*compartment=\")(.*?)(\"\\s.*)", "$2");
					Compound c = new Compound(id, name, compartment);
					ps.addPrecursor(c);
				}
				if(line.startsWith("</precursorSet>")){ // new precursorSet
					precursorSetSolutions.add(ps);
					ps = new PrecursorSet(0, 0, 0);
				}
				line = br.readLine();
			}
			br.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return precursorSetSolutions;
	}
}
