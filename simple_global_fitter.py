import os, sys
from galaxy import eggs
from galaxy.tools.util import hyphy_util



#Maybe it's better to put this into galaxy.tools.util.hyphy_util like similar stuff
def get_simple_global_fitter_config_filename(input_filename, output_filename, code, model, bf_file):
    contents = """
_genomeScreenOptions = {};

/* all paths are either absolute or relative to the DATA READER;
   the first 4 options are for the data reader*/

_genomeScreenOptions ["00"] = "../AnalysisModules/SimpleGlobalFitter.bf";
/* which analysis to run on each gene; */
_genomeScreenOptions ["01"] = "../Writers/TAB.bf";
/* what output to produce; */
_genomeScreenOptions ["02"] = "";
/* tree string; leave blank to use NJ for inference in some analyses */
/* try to paste "(human, chimp, mouse)" to see the use of a fixed tree */
_genomeScreenOptions ["03"] = "%s";
/* alignment file */
_genomeScreenOptions ["04"] = "%s";
/* output csv file */

/* options for the analysis */
_genomeScreenOptions ["05"] = "%s";
/* genetic code */
_genomeScreenOptions ["06"] = "%s";
/* nucleotide bias string; can define any of the 203 models */

ExecuteAFile ("%s", _genomeScreenOptions);
""" % (input_filename, output_filename, code, model, bf_file)
    return hyphy_util.get_filled_temp_filename(contents)



#Retrieve hyphy path, this will need to be the same across the cluster
tool_data = sys.argv.pop()
HYPHY_PATH = os.path.join(tool_data, "HYPHY")
HYPHY_EXECUTABLE = os.path.join(HYPHY_PATH, "HYPHY")

#get config file
bf_file = os.path.dirname(__file__) + "/GenomeFitters/DataReaders/FastaReader.bf"

#Read command line arguments
input_filename = os.path.abspath(sys.argv[1].strip())
output_filename = os.path.abspath(sys.argv[2].strip())
code = sys.argv[3].strip()
model = sys.argv[4].strip()

#setup Config file
config_filename = get_simple_global_fitter_config_filename(input_filename, output_filename, code, model, bf_file)

#Run Hyphy
hyphy_cmd = "%s BASEPATH=%s USEPATH=/dev/null %s" % (HYPHY_EXECUTABLE, HYPHY_PATH, config_filename)
hyphy = os.popen(hyphy_cmd, 'r')
#print hyphy.read()
hyphy.close()

#remove temporary files
os.unlink(config_filename)

