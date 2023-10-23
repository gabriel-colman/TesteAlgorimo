package utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import metabolicNetwork.Compound;
import metabolicNetwork.MetabolicNetwork;
import metabolicNetwork.Reaction;

public class MetabolicNetworkSBMLWriter 
{
	String filename;
	File networkFile;
	BufferedWriter bw;
	
	public MetabolicNetworkSBMLWriter(String filename) 
	{
		this.filename = filename;
	}
	
	private void closeFile() {
		if (networkFile != null) {
			try {
				bw.close();
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(1);
			}
			bw = null;
			networkFile = null;
		}
	}
	
	public void write(MetabolicNetwork network) 
	{
		networkFile = new File(filename);
		try {
			bw = new BufferedWriter(new FileWriter(networkFile));
			writeHeader();
			writeListOfCompounds(network);
			writeListOfReactions(network);
			writeTail();
			closeFile();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}	
	}

	private void writeHeader() throws IOException {
		bw.write("<?xml version=\"1.0\"  encoding=\"UTF-8\"?>");
		bw.write("\n<sbml xmlns=\"http://www.sbml.org/sbml/level2\" version=\"1\" level=\"2\" xmlns:html=\"http://www.w3.org/1999/xhtml\">");
		bw.write("\n<model id=\"automatically_generated\" name=\"auto\">");
	}
	
	private void writeListOfCompounds(MetabolicNetwork network) throws IOException {
		bw.write("\n<listOfSpecies>");
//		  <species id="M_abt_e" name="M_L_Arabinitol_C5H12O5" initialAmount="0" compartment="cytoplasm" boundaryCondition="false" />
		for(Compound c: network.getCompounds().values()) {
			String name = "";
			String id = "";
			if (c.getId() != null){
				id = c.getId();
			}
			if (c.getName() != null){
				name = c.getName();
			}
			bw.write("\n\t<species id=\""+StringUtils.sbmlEncode(id)+"\" name=\""+StringUtils.sbmlEncode(name)+"\""+"/>");
		}
		bw.write("\n</listOfSpecies>");	
	}
	
	private void writeListOfReactions(MetabolicNetwork network) throws IOException {
		bw.write("\n<listOfReactions>");
//		  <reaction id="R_PANTS" name="R_pantothenate_synthase" reversible="false">
//		      <listOfReactants>
//		        <speciesReference species="M_pant_R_c" stoichiometry="1.000000"/>
//		        <speciesReference species="M_ala_B_c" stoichiometry="1.000000"/>
//		      </listOfReactants>
//		      <listOfProducts>
//		        <speciesReference species="M_pnto_R_c" stoichiometry="1.000000"/>
//		      </listOfProducts>
//		    </reaction>
		for(Reaction r: network.getReactions().values()) {

			String rName = (r.getName()==null)?StringUtils.sbmlEncode(r.getId()):StringUtils.sbmlEncode(r.getName());
			bw.write("\n\t<reaction id=\""+StringUtils.sbmlEncode(r.getId())+"\" name=\""+rName+"\" reversible=\""+r.isReversible()+"\">");

			bw.write("\n\t\t<listOfReactants>");
			for(Compound sub: r.getSubstrates().values()) {
				bw.write("\n\t\t\t<speciesReference species=\""+StringUtils.sbmlEncode(sub.getId())+"\" stoichiometry=\"" + (double) r.getSubstrateStochiometricValue(sub) + "\"/>");				
			}
			bw.write("\n\t\t</listOfReactants>");

			bw.write("\n\t\t<listOfProducts>");
			for(Compound prod: r.getProduces().values()) {
				bw.write("\n\t\t\t<speciesReference species=\""+StringUtils.sbmlEncode(prod.getId())+"\" stoichiometry=\"" + (double) r.getProductStochiometricValue(prod) + "\"/>");				
			}
			bw.write("\n\t\t</listOfProducts>");

			bw.write("\n\t</reaction>");			
		}
		bw.write("\n</listOfReactions>");			
	}
	
	private void writeTail() throws IOException {
		bw.write("\n</model>");
		bw.write("\n</sbml>");
	}
    
}
