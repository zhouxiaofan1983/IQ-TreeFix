#
# Python module for SH test using IQTREE
# Modified from the original RAxML model: raxmlmodel.py 
# by Xiaofan Zhou, 2015/09/18
#

# treefix libraries
from treefix.models import StatModel

# python libraries
import os, sys, re
import optparse
import tempfile

# rasmus libraries
from rasmus import util

# compbio libraries
from compbio import phylip

class CoarseModelSpeed(StatModel):
    """Computes test statistics using IQTREE"""

    def __init__(self, extra):
        """Initializes the IQTREE model"""
        StatModel.__init__(self, extra)

        self.VERSION = "1.0"
		
        parser = optparse.OptionParser(prog="IQTREEModel")
        parser.add_option("-c", "--cpu", dest="cpu",
                          metavar="<cpu>",
                          default=1,
                          help="number of threads to use (default: 1)")
        parser.add_option("-m", "--model", dest="model",
                          metavar="<model>",
                          default="GTR+G",
                          help="model of DNA, AA, or CODON substitution (default: GTR+G)")
        parser.add_option("-t", "--type", dest="type",
                          metavar="<type>",
                          default="DNA",
                          help="data type, one of the following: DNA, AA, NT2AA, CODON (default: DNA)")
        parser.add_option("-r", "--replicates", dest="rep",
                          metavar="<rep>",
                          default=10000,
                          help="number of RELL replicates for the tree topology test (default: 10000)")
        parser.add_option("-p", "--prefix", dest="pre",
                          metavar="<pre>",
                          default="",
                          help="the prefix for temporary output files (default: '')")
        self.parser = parser

        StatModel._parse_args(self, extra)

    def __del__(self):
        """Cleans up the IQTREE model"""
        
        os.remove(self.btreefile)
        os.remove(self.seqfile)
        if self.type == "DNA":
            os.remove(self.modelfile)

    def optimize_model(self, gtree, aln):
        """Optimizes the IQTREE model"""
        
	fd, btreefile = tempfile.mkstemp('.btree')
        os.close(fd)
	gtree.write(btreefile)
        self.btreefile = btreefile

        fd, seqfile = tempfile.mkstemp('.align')
        os.close(fd)
        out = util.open_stream(seqfile, "w")
        phylip.write_phylip_align(out, aln, strip_names=False)
        out.close()
        self.seqfile = seqfile
	
        # optimize all model parameters and create a model file for DNA data, or optimize the gamma shape parameter for other data type if necessary
	if self.type == "DNA" or "+G" in self.model:
            # optimize model parameters on user specified tree
            os.system('iqtree-omp -nt %s -m %s -st %s -s %s -te %s -pre %s.treefix_tmp > /dev/null' % (self.cpu, self.model, self.type, self.seqfile, self.btreefile, self.pre))
            # DNA data
            if self.type == "DNA":
                # create model file for DNA data
                fd, modelfile = tempfile.mkstemp('.model')
                os.close(fd)
                self.modelfile = modelfile
                # read estimated model parameters from IQTREE output
                fi = open("%s.treefix_tmp.iqtree" % self.pre, 'r')
                fo = open(modelfile, 'w')
                line = fi.readline()
                while line:
                    if "Rate parameter R:" in line:
                        fi.readline()
                        while True:
                            line = fi.readline().rstrip()
                            if ":" in line or "=" in line:
                                param = re.split('\s+', line)
                                fo.write("%s\n" % param[-1])
                            else:
                                break
                    elif "State frequencies:" in line:
                        if "equal frequencies" in line:
                            for i in range(4):
                                fo.write("0.25\n")
                        else:
                            fi.readline()
                            while True:
                                line = fi.readline().rstrip()
                                if ":" in line or "=" in line:
                                    param = re.split('\s+', line)
                                    fo.write("%s\n" % param[-1])
                                else:
                                    break
                    elif "Gamma shape alpha:" in line:
                        param = re.split('\s+', line)
                        self.gamma = param[3]
                        break
                    line = fi.readline()
                fi.close()
                fo.close()
                # set command for stats
                if "+G" in self.model:
                    self.cmd = 'iqtree-omp -nt %s -m %s+G+Fu -a %s -st %s -s %s -te %s -z %s.treefix_tmp.tree -zb %s -pre %s.treefix_tmp > /dev/null' % (self.cpu, self.modelfile, self.gamma, self.type, self.seqfile, self.btreefile, self.pre, self.rep, self.pre)
                else:
                    self.cmd = 'iqtree-omp -nt %s -m %s+Fu -st %s -s %s -te %s -z %s.treefix_tmp.tree -zb %s -pre %s.treefix_tmp > /dev/null' % (self.cpu, self.modelfile, self.type, self.seqfile, self.btreefile, self.pre, self.rep, self.pre)
            # other data type
            elif "+G" in self.model:
                # read estimated model parameters from IQTREE outptut
                fi = open("%s.treefix_tmp.iqtree" % self.pre, 'r')
                for line in fi:
                    if "Gamma shape alpha:" in line:
                        param = re.split('\s+', line)
                        self.gamma = param[3]
                        break
                fi.close()
                # set command for stats
                self.cmd = 'iqtree-omp -nt %s -m %s -a %s -st %s -s %s -te %s -z %s.treefix_tmp.tree -zb %s -pre %s.treefix_tmp > /dev/null' % (self.cpu, self.model, self.gamma, self.type, self.seqfile, self.btreefile, self.pre, self.rep, self.pre)
            # clean up files generated during model optimization
            os.system('rm %s.treefix_tmp.*' % self.pre)
        else:
            # set command for stats
            self.cmd = 'iqtree-omp -nt %s -m %s -st %s -s %s -te %s -z %s.treefix_tmp.tree -zb %s -pre %s.treefix_tmp > /dev/null' % (self.cpu, self.model, self.type, self.seqfile, self.btreefile, self.pre, self.rep, self.pre)

    def compute_lik_test(self, gtree, stat="SH", alternative=None):
        """Computes the test statistic 'stat' using IQTREE"""
	
	pval = 0
	Dlnl = 0
        
	fd, gtreefile = tempfile.mkstemp('.gtree')
        os.close(fd)
        gtree.write(gtreefile)
	os.system('cat %s %s > %s.treefix_tmp.tree' % (self.btreefile, gtreefile, self.pre))
        
        os.system('%s' % self.cmd)
        
        f = open('%s.treefix_tmp.iqtree' % self.pre, 'r')
        line = f.readline()
	while line:
            if "USER TREES" in line:
                for i in range(8):
	            line = f.readline()
                stats = re.split('\s+', line)
                pval = float(stats[8])
                Dlnl = float(stats[3])
                break
            line = f.readline()
        f.close()        

        os.remove(gtreefile)
	os.system('rm %s.treefix_tmp.*' % self.pre)
        
	return pval, Dlnl

class CoarseModel(StatModel):
    """Computes test statistics using IQTREE"""

    def __init__(self, extra):
        """Initializes the IQTREE model"""
        StatModel.__init__(self, extra)

        self.VERSION = "1.0"
		
        parser = optparse.OptionParser(prog="IQTREEModel")
        parser.add_option("-c", "--cpu", dest="cpu",
                          metavar="<cpu>",
                          default=1,
                          help="number of threads to use (default: 1)")
        parser.add_option("-m", "--model", dest="model",
                          metavar="<model>",
                          default="GTR+G",
                          help="model of DNA, AA, or CODON substitution (default: GTR+G)")
        parser.add_option("-t", "--type", dest="type",
                          metavar="<type>",
                          default="DNA",
                          help="data type, one of the following: DNA, AA, NT2AA, CODON (default: DNA)")
        parser.add_option("-r", "--replicates", dest="rep",
                          metavar="<rep>",
                          default=10000,
                          help="number of RELL replicates for the tree topology test (default: 10000)")
        parser.add_option("-p", "--prefix", dest="pre",
                          metavar="<pre>",
                          default="",
                          help="the prefix for temporary output files (default: '')")
        self.parser = parser

        StatModel._parse_args(self, extra)

    def __del__(self):
        """Cleans up the IQTREE model"""
        
        os.remove(self.btreefile)
        os.remove(self.seqfile)

    def optimize_model(self, gtree, aln):
        """Optimizes the IQTREE model"""
         
	fd, btreefile = tempfile.mkstemp('.btree')
        os.close(fd)
	gtree.write(btreefile)
        self.btreefile = btreefile
        
        fd, seqfile = tempfile.mkstemp('.align')
        os.close(fd)
        out = util.open_stream(seqfile, "w")
        phylip.write_phylip_align(out, aln, strip_names=False)
        out.close()
        self.seqfile = seqfile

    def compute_lik_test(self, gtree, stat="SH", alternative=None):
        """Computes the test statistic 'stat' using IQTREE"""
	
	pval = 0
	Dlnl = 0
	
	fd, gtreefile = tempfile.mkstemp('.gtree')
        os.close(fd)
        gtree.write(gtreefile)
	os.system('cat %s %s > %s.treefix_tmp.tree' % (self.btreefile, gtreefile, self.pre))
       
	os.system('iqtree-omp -nt %s -m %s -st %s -s %s -te %s -z %s.treefix_tmp.tree -zb %s -pre %s.treefix_tmp > /dev/null' % (self.cpu, self.model, self.type, self.seqfile, self.btreefile, self.pre, self.rep, self.pre))

        f = open('%s.treefix_tmp.iqtree' % self.pre, 'r')
        line = f.readline()
	while line:
            if "USER TREES" in line:
                for i in range(8):
	            line = f.readline()
                stats = re.split('\s+', line)
                pval = float(stats[8])
                Dlnl = float(stats[3])
                break
            line = f.readline()
        f.close()        
        
        os.remove(gtreefile)
	os.system('rm %s.treefix_tmp.*' % self.pre)

	return pval, Dlnl

class FineModelSpeed(StatModel):
    """Computes test statistics using IQTREE"""

    def __init__(self, extra):
        """Initializes the IQTREE model"""
        StatModel.__init__(self, extra)

        self.VERSION = "1.0"
		
        parser = optparse.OptionParser(prog="IQTREEModel")
        parser.add_option("-c", "--cpu", dest="cpu",
                          metavar="<cpu>",
                          default=1,
                          help="number of threads to use (default: 1)")
        parser.add_option("-m", "--model", dest="model",
                          metavar="<model>",
                          default="GTR+G",
                          help="model of DNA, AA, or CODON substitution (default: GTR+G)")
        parser.add_option("-t", "--type", dest="type",
                          metavar="<type>",
                          default="DNA",
                          help="data type, one of the following: DNA, AA, NT2AA, CODON (default: DNA)")
        parser.add_option("-r", "--replicates", dest="rep",
                          metavar="<rep>",
                          default=1,
                          help="number of RELL replicates for the tree topology test (default: 1)")
        parser.add_option("-p", "--prefix", dest="pre",
                          metavar="<pre>",
                          default="",
                          help="the prefix for temporary output files (default: '')")
        self.parser = parser

        StatModel._parse_args(self, extra)
    
    def __del__(self):
        """Cleans up the IQTREE model"""
        
        os.remove(self.seqfile)
        if self.type == "DNA":
            os.remove(self.modelfile)

    def optimize_model(self, gtree, aln):
        """Optimizes the IQTREE model"""
	
	fd, btreefile = tempfile.mkstemp('.btree')
        os.close(fd)
	gtree.write(btreefile)

        fd, seqfile = tempfile.mkstemp('.align')
        os.close(fd)
        out = util.open_stream(seqfile, "w")
        phylip.write_phylip_align(out, aln, strip_names=False)
        out.close()
	self.seqfile = seqfile

	fd, bsitelhfile = tempfile.mkstemp('.bsitelh')
        os.close(fd)
	
        # optimize model parameters on user specified tree
        os.system('iqtree-omp -nt %s -m %s -st %s -s %s -te %s -pre %s.treefix_tmp -wsl > /dev/null' % (self.cpu, self.model, self.type, self.seqfile, btreefile, self.pre))
        # optimize all model parameters and create a model file for DNA data, or optimize the gamma shape parameter for other data type if necessary
	if self.type == "DNA" or "+G" in self.model:
            # DNA data
            if self.type == "DNA":
                # create model file for DNA data
                fd, modelfile = tempfile.mkstemp('.model')
                os.close(fd)
                self.modelfile = modelfile
                # read estimated model parameters from IQTREE output
                fi = open("%s.treefix_tmp.iqtree" % self.pre, 'r')
                fo = open(modelfile, 'w')
                line = fi.readline()
                while line:
                    if "Rate parameter R:" in line:
                        fi.readline()
                        while True:
                            line = fi.readline().rstrip()
                            if ":" in line or "=" in line:
                                param = re.split('\s+', line)
                                fo.write("%s\n" % param[-1])
                            else:
                                break
                    elif "State frequencies:" in line:
                        if "equal frequencies" in line:
                            for i in range(4):
                                fo.write("0.25\n")
                        else:
                            fi.readline()
                            while True:
                                line = fi.readline().rstrip()
                                if ":" in line or "=" in line:
                                    param = re.split('\s+', line)
                                    fo.write("%s\n" % param[-1])
                                else:
                                    break
                    elif "Gamma shape alpha:" in line:
                        param = re.split('\s+', line)
                        self.gamma = param[3]
                        break
                    line = fi.readline()
                fi.close()
                fo.close()
                # set command for stats
                if "+G" in self.model:
                    self.cmd = 'iqtree-omp -nt %s -m %s+G+Fu -a %s -st %s -s %s' % (self.cpu, self.modelfile, self.gamma, self.type, self.seqfile)
                else:
                    self.cmd = 'iqtree-omp -nt %s -m %s+Fu -st %s -s %s' % (self.cpu, self.modelfile, self.type, self.seqfile)
            # other data type
            elif "+G" in self.model:
                # read estimated model parameters from IQTREE outptut
                fi = open("%s.treefix_tmp.iqtree" % self.pre, 'r')
                for line in fi:
                    if "Gamma shape alpha:" in line:
                        param = re.split('\s+', line)
                        self.gamma = param[3]
                        break
                fi.close()
                # set command for stats
                self.cmd = 'iqtree-omp -nt %s -m %s -a %s -st %s -s %s' % (self.cpu, self.model, self.gamma, self.type, self.seqfile)
        else:
            # set command for stats
            self.cmd = 'iqtree-omp -nt %s -m %s -st %s -s %s' % (self.cpu, self.model, self.type, self.seqfile)
        
        f = open("%s.treefix_tmp.sitelh" % self.pre, 'r')
        self.bsitelh = f.readline().replace("1", "2", 1) + f.readline().replace("Site_Lh", "Tree1", 1)
        f.close()

	os.system('rm %s.treefix_tmp.*' % self.pre)
        os.remove(btreefile)

    def compute_lik_test(self, gtree, stat="SH", alternative=None):
        """Computes the test statistic 'stat' using IQTREE"""
	
	pval = 0
	Dlnl = 0
 
	fd, gtreefile = tempfile.mkstemp('.gtree')
        os.close(fd)
        gtree.write(gtreefile)
      
	os.system('%s -te %s -pre %s.treefix_tmp -wsl > /dev/null' % (self.cmd, gtreefile, self.pre))
        fo = open("%s.sitelh" % self.pre, 'w')
        fo.write(self.bsitelh)
	fi = open("%s.treefix_tmp.sitelh" % self.pre, 'r')
        fi.readline()
        fo.write(fi.readline().replace("Site_Lh", "Tree2", 1))
        fi.close()
        fo.close()
	os.system('seqmt --puzzle %s.sitelh > /dev/null' % self.pre)
	if stat == "AU":
	    os.system('makermt -b %s %s > /dev/null' % (self.rep, self.pre))
	else:
	    os.system('makermt -f -b %s %s > /dev/null' % (self.rep, self.pre))
	os.system('consel %s > /dev/null' % self.pre)
        os.system('catpv -s 1 %s.pv > %s.catpv.out' % (self.pre, self.pre))
 
	f = open("%s.catpv.out" % self.pre, 'r')
        for i in range(6):
            line = f.readline()
        stats = re.split('\s+', line)
        if stat == "AU":
            pval = float(stats[4])
        else:
            pval = float(stats[10])
        Dlnl = float(stats[3])

        os.remove(gtreefile)
	os.system('rm %s.treefix_tmp.* %s.sitelh %s.mt %s.rmt %s.vt %s.pv %s.catpv.out' % (self.pre, self.pre, self.pre, self.pre, self.pre, self.pre, self.pre))

	return pval, Dlnl

class FineModel(StatModel):
    """Computes test statistics using IQTREE"""

    def __init__(self, extra):
        """Initializes the IQTREE model"""
        StatModel.__init__(self, extra)

        self.VERSION = "1.0"
		
        parser = optparse.OptionParser(prog="IQTREEModel")
        parser.add_option("-c", "--cpu", dest="cpu",
                          metavar="<cpu>",
                          default=1,
                          help="number of threads to use (default: 1)")
        parser.add_option("-m", "--model", dest="model",
                          metavar="<model>",
                          default="GTR+G",
                          help="model of DNA, AA, or CODON substitution (default: GTR+G)")
        parser.add_option("-t", "--type", dest="type",
                          metavar="<type>",
                          default="DNA",
                          help="data type, one of the following: DNA, AA, NT2AA, CODON (default: DNA)")
        parser.add_option("-r", "--replicates", dest="rep",
                          metavar="<rep>",
                          default=1,
                          help="number of RELL replicates for the tree topology test (default: 1)")
        parser.add_option("-p", "--prefix", dest="pre",
                          metavar="<pre>",
                          default="",
                          help="the prefix for temporary output files (default: '')")
        self.parser = parser

        StatModel._parse_args(self, extra)
    
    def __del__(self):
        """Cleans up the IQTREE model"""
        
        os.remove(self.seqfile)

    def optimize_model(self, gtree, aln):
        """Optimizes the IQTREE model"""
	
	fd, btreefile = tempfile.mkstemp('.btree')
        os.close(fd)
	gtree.write(btreefile)

        fd, seqfile = tempfile.mkstemp('.align')
        os.close(fd)
        out = util.open_stream(seqfile, "w")
        phylip.write_phylip_align(out, aln, strip_names=False)
        out.close()
	self.seqfile = seqfile

	fd, bsitelhfile = tempfile.mkstemp('.bsitelh')
        os.close(fd)
	
	os.system('iqtree-omp -nt %s -m %s -st %s -s %s -te %s -pre %s.treefix_tmp -wsl > /dev/null' % (self.cpu, self.model, self.type, self.seqfile, btreefile, self.pre))
        
        f = open("%s.treefix_tmp.sitelh" % self.pre, 'r')
        self.bsitelh = f.readline().replace("1", "2", 1) + f.readline().replace("Site_Lh", "Tree1", 1)
        f.close()

	os.system('rm %s.treefix_tmp.*' % self.pre)
        os.remove(btreefile)

    def compute_lik_test(self, gtree, stat="SH", alternative=None):
        """Computes the test statistic 'stat' using IQTREE"""
	
	pval = 0
	Dlnl = 0
 
	fd, gtreefile = tempfile.mkstemp('.gtree')
        os.close(fd)
        gtree.write(gtreefile)
      
	os.system('iqtree-omp -nt %s -m %s -st %s -s %s -te %s -pre %s.treefix_tmp -wsl > /dev/null' % (self.cpu, self.model, self.type, self.seqfile, gtreefile, self.pre))
        fo = open("%s.sitelh" % self.pre, 'w')
        fo.write(self.bsitelh)
	fi = open("%s.treefix_tmp.sitelh" % self.pre, 'r')
        fi.readline()
        fo.write(fi.readline().replace("Site_Lh", "Tree2", 1))
        fi.close()
        fo.close()
	os.system('seqmt --puzzle %s.sitelh > /dev/null' % self.pre)
	if stat == "AU":
	    os.system('makermt -b %s %s > /dev/null' % (self.rep, self.pre))
	else:
	    os.system('makermt -f -b %s %s > /dev/null' % (self.rep, self.pre))
	os.system('consel %s > /dev/null' % self.pre)
        os.system('catpv -s 1 %s.pv > %s.catpv.out' % (self.pre, self.pre))
 
	f = open("%s.catpv.out" % self.pre, 'r')
        for i in range(6):
            line = f.readline()
        stats = re.split('\s+', line)
        if stat == "AU":
            pval = float(stats[4])
        else:
            pval = float(stats[10])
        Dlnl = float(stats[3])

        os.remove(gtreefile)
	os.system('rm %s.treefix_tmp.* %s.sitelh %s.mt %s.rmt %s.vt %s.pv %s.catpv.out' % (self.pre, self.pre, self.pre, self.pre, self.pre, self.pre, self.pre))

	return pval, Dlnl
