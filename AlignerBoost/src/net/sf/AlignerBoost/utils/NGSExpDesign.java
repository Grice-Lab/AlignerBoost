/**
 * A helper class to parse tab-delimited NGS experimental design files
 */
package net.sf.AlignerBoost.utils;
import java.io.*;
import java.util.*;
import java.util.regex.*;

/**
 * @author Qi Zheng
 * @version 1.1.1
 * @since 1.1
 */
public class NGSExpDesign {
	/**
	 * Construct a NGSExpDesign instance from a given file
	 * @param designFileName experimental design filename
	 * @throws IOException throws IOException if paramFileName does not exist or other IO errors happened
	 */
	public NGSExpDesign(String designFileName) throws IOException, IllegalArgumentException {
		// initiate fields
		globalOpts = new HashMap<String, String>();
		libNames = new ArrayList<String>();
		optNames = new ArrayList<String>();
		libOpts = new HashMap<String, Map<String, String>>();
		
		// read expDesign file
		BufferedReader in = new BufferedReader(new FileReader(designFileName));
		String line = null;
		while((line = in.readLine()) != null) {
			if(line.startsWith("#")) { // a header line
				Matcher match = globalPat.matcher(line);
				if(match.find()) // a global opt line
					globalOpts.put(match.group(1),  match.group(2));
				match = localPat.matcher(line);
				if(match.find()) // a local opt line
					optNames.add(match.group(1));
			}
			else { // a opt-value line
				String[] optValues = line.split("\t");
				if(optNames.size() < optValues.length)
					throw new IllegalArgumentException("Incorrect number of value fields found at\n$line\nNo more than " +
							optNames.size() + " fields are allowed");
				String libName = optValues[0];
				if(libOpts.containsKey(libName))
					throw new IllegalArgmentException("Redundant libname found for '" + libName + "'");
				else
					libOpts.put(libName, new HashMap<String, String>()); // Init inner Map as a HashMap
				// pair opt names with values
				for(int i = 0; i < optNames.size(); i++) {
					String val = i < optValues.length ? optValues[i] : "";
					libOpts.get(libName).put(optNames.get(i), optValues[i]);
				}
				// fill default value for this libOpts
				setDefaultOptValues(libOpts.get(libName));
			}	
		}
		
	}
	
	private void setDefaultOptValues(Map<String, String> opts) {
		if(opts.get("read_file").equals("NA"))
			opts.put("read_file", opts.get("libname") + ".fastq");
		if(opts.get("ref_index").equals("NA"))
			opts.put("ref_index", opts.get("ref_genome"));
		if(opts.get("transcriptome_index").equals("NA"))
			opts.put("transcriptome_index", "transcriptome/" + opts.get("transcriptome_GFF").replaceFirst("(?i:\\.gff)", ""));		
			
		/*				$lib_opt{$libname}{'ref_index'} = $lib_opt{$libname}{'ref_genome'} if($lib_opt{$libname}{'ref_index'} eq 'NA');
		if($lib_opt{$libname}{'transcriptome_index'} eq 'NA') {
			$lib_opt{$libname}{'transcriptome_index'} = 'transcriptome/' . $lib_opt{$libname}{'transcriptome_GFF'};
			$lib_opt{$libname}{'transcriptome_index'} =~ s/\.gff$//i;
		}
		if(!$lib_opt{$libname}{'do_trim'} && !$lib_opt{$libname}{'has_spliced'}) { # fix min_insert as read_len
			$lib_opt{$libname}{'min_insert'} = $lib_opt{$libname}{'read_len'};
		}*/
	}

	private Map<String, String> globalOpts;
	private List<String> libNames;
	private List<String> optNames;
	private Map<String, Map<String, String>> libOpts;
	
	private Pattern globalPat = Pattern.compile("^## (\\w+)=(.*):");
	private Pattern localPat = Pattern.compile("^### (\\w+):");
}
