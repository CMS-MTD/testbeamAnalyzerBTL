# /usr/bin/python

#Author: Ben Tannenwald
#Date: April 21, 2019
#Purpose: Script to submit condor jobs for april TB timing analysis

import os,sys, argparse, getpass

# *** 0. setup parser for command line
parser = argparse.ArgumentParser()
parser.add_argument("--bar", help="which bar to process, e.g. box1, box2, box3, single")
parser.add_argument("--firstRun", help="first run to analyze")
parser.add_argument("--lastRun", help="last run to analyze")
parser.add_argument("--biasVoltage", help="bias voltage [V]")
parser.add_argument("--timeAlgo", help="time fitting algorithm", default="IL_50")
args = parser.parse_args()

if len(vars(args)) != 5 or len(sys.argv)==1: # 4 --> four: one for each options
    os.system('python analyzeBarForTiming.py -h')
    quit()

# ** A. Test bar option and exit if unexpected
if(args.bar is None):
    print "#### Need to bar to analyze using --bar <box1/box2/box3/single> ####\nEXITING"
    quit()
else:
    if (args.bar == "box1" or args.bar == "box2" or args.bar == "box3" or args.bar == "single")==False:
        print "#### Need to bar to analyze using --bar <box1/box2/box3/single> ####\nEXITING"
        quit()
    else:
        print '-- Setting  bar = {0}'.format(args.bar)

# ** B. first run
if(args.firstRun is None):
    print "#### Need to specify first run to analyze --firstRun <run number###\nEXITING"
    quit()
else:
    if (args.firstRun).isdigit == False:
        print "#### Need to specify first run to analyze --firstRun <run number###\nEXITING"
        quit()
    else:
        print '-- Setting firstRun = {0}'.format(args.firstRun)

# ** C. last run
if(args.lastRun is None):
    print "#### Need to specify last run to analyze --lastRun <run number###\nEXITING"
    quit()
else:
    if (args.lastRun).isdigit == False:
        print "#### Need to specify last run to analyze --lastRun <run number###\nEXITING"
        quit()
    else:
        print '-- Setting lastRun = {0}'.format(args.lastRun)

# ** D. time Algolast run
if(args.timeAlgo == "IL_50"):
    print '-- Default timeAlgo = {0}'.format(args.timeAlgo)
else:
    print '-- Setting lastRun = {0}'.format(args.timeAlgo)

# ** E. bias voltage
if(args.biasVoltage is None):
    print "#### Need to specify bias voltage to analyze --biasVoltage <bias voltage> ###\nEXITING"
    quit()
else:
    if (args.lastRun).isdigit == False:
        print "#### Need to specify bias voltage to analyze --biasVoltage <bias voltage> ###\nEXITING"
        quit()
    else:
        print '-- Setting biasVoltage = {0}'.format(args.biasVoltage)


# ** F. Exit if no grid proxy
if ( not os.path.exists(os.path.expandvars("$X509_USER_PROXY")) ):
    print "#### No GRID PROXY detected. Please do voms-proxy-init -voms cms before submitting Condor jobs ####.\nEXITING"
    quit()



# *** 1. Parse together output directory and create if does not exist
# ** A. Super-Top level directory for testbeam analysis
username = getpass.getuser()
eosDir = "/eos/uscms/store/user/{0}/testbeam_04-2019/".format(username)
if ( not os.path.isdir(eosDir) ):
    print "Specified top directory {0} DNE.\nCREATING NOW".format(eosDir)
    os.system("rm -rf {0}".format(outputDir))
    os.system("mkdir {0}".format(eosDir))

# ** B. Top level directory of which runs
outputDir = 'Runs_'+args.firstRun+'to'+args.lastRun
if ( not os.path.isdir(outputDir) ):
    print "Specified run directory {0} DNE.\nCREATING NOW".format(outputDir)
    os.system("rm -rf {0}".format(outputDir))
    os.system("mkdir {0}".format(outputDir))
if ( not os.path.isdir(eosDir+'/{0}/'.format(outputDir)) ):
    os.system("mkdir "+eosDir+'/{0}/'.format(outputDir))

# ** C. Sub-level directory of what bias
outputDir = outputDir + '/' + args.biasVoltage + 'V'
if ( not os.path.isdir(outputDir) ):
    print "Specified bias voltage sub-directory {0} DNE.\nCREATING NOW".format(outputDir)
    os.system("rm -rf {0}".format(outputDir))
    os.system("mkdir {0}".format(outputDir))
if ( not os.path.isdir(eosDir+'/{0}/'.format(outputDir)) ):
    os.system("mkdir "+eosDir+'/{0}/'.format(outputDir))

# ** D. Sub-level directory of which bar
outputDir = outputDir + '/' + args.bar
if ( not os.path.isdir(outputDir) ):
    print "Specified bar sub-directory {0} DNE.\nCREATING NOW".format(outputDir)
    os.system("rm -rf {0}".format(outputDir))
    os.system("mkdir {0}".format(outputDir))
if ( not os.path.isdir(eosDir+'/{0}/'.format(outputDir)) ):
    os.system("rm -rf "+eosDir+'/{0}/'.format(outputDir))
    os.system("mkdir "+eosDir+'/{0}/'.format(outputDir))

# ** D. Sub-level directory of which time algorithm
outputDir = outputDir + '/' + args.timeAlgo 
if ( not os.path.isdir(outputDir) ):
    print "Specified timeAlgo sub-directory {0} DNE.\nCREATING NOW".format(outputDir)
    os.system("rm -rf {0}".format(outputDir))
    os.system("mkdir {0}".format(outputDir))
if ( not os.path.isdir(eosDir+'/{0}/'.format(outputDir)) ):
    print "EOS Specified timeAlgo sub-directory {0} DNE.\nCREATING NOW".format(outputDir)
    os.system("rm -rf "+eosDir+'/{0}'.format(outputDir))
    os.system("mkdir "+eosDir+'/{0}/'.format(outputDir))

# ** E. Make folders for condor output storage
if ( not os.path.isdir( (outputDir + '/condor_logs/') ) ):
    os.system("mkdir {0}".format( (outputDir + '/condor_logs/')) )
if ( not os.path.isdir( (outputDir + '/condor_err/') ) ):
    os.system("mkdir {0}".format( (outputDir + '/condor_err/')) )
if ( not os.path.isdir( (outputDir + '/condor_out/') ) ):
    os.system("mkdir {0}".format( (outputDir + '/condor_out/')) )
    
print '-- Setting outputDir = {0}'.format(outputDir)




# *** 2. Create .tar of directory and store in personal EOS
print "##########     Tarring workdir     ##########"
tarball_name = "{0}.tar.gz".format(outputDir.replace('/', '__'))
# remove EOS tarball
if os.path.isfile( eosDir+'/'+outputDir+'/'+tarball_name ) :
    os.system('rm {0}'.format(eosDir+'/'+outputDir+'/'+tarball_name))
    print ("WARNING: Removing {0}.\nTHIS MEANS YOU HAVE RUN THIS CONFIGURATION ONCE BEFORE. IF YOU HAVEN'T SAVED THE PREVIOUS TARBALL BEFORE, IT'S TOO LATE NOW.\n".format(eosDir+'/'+outputDir+'/'+tarball_name))
# remove local tarball
if os.path.isfile( './'+tarball_name ) :
    os.system('rm {0}'.format(tarball_name)) 
os.system("tar -cvzf {0} ./ --exclude 'Runs*' --exclude 'submitOneFile_' --exclude '*.tar.gz' --exclude '*.*~' --exclude 'temp' --exclude 'test*' --exclude '*04-*-19' --exclude '.git'".format(tarball_name))
os.system("xrdcp {0} root://cmseos.fnal.gov//store/user/{1}/testbeam_04-2019/{2}/".format(tarball_name, username, outputDir))



# *** 3. Create temporary .pdl file for condor submission
print "\n##########     Submitting Condor jobs     ##########\n"
jdl_filename = "submitOneFile_{0}.jdl".format(outputDir.replace('/', '__'))

os.system("touch {0}".format(jdl_filename))
os.system("echo universe = vanilla > {0}".format(jdl_filename))
os.system("echo Executable = runOneFile.csh >> {0}".format(jdl_filename))
os.system("echo Should_Transfer_Files = YES >> {0}".format(jdl_filename))
os.system("echo WhenToTransferOutput = ON_EXIT >> {0}".format(jdl_filename))
os.system("echo Transfer_Input_Files = runOneFile.csh, {0} >> {1}".format(tarball_name, jdl_filename))
os.system("echo Output = {0}/condor_out/outfile_{1}.out  >> {2}".format(outputDir, outputDir.replace('/', '__'), jdl_filename))
os.system("echo Error = {0}/condor_err/outfile_{1}.err >> {2}".format(outputDir, outputDir.replace('/', '__'), jdl_filename))
os.system("echo Log = {0}/condor_logs/outfile_{1}.log >> {2}".format(outputDir, outputDir.replace('/', '__'), jdl_filename))
os.system("echo x509userproxy = ${{X509_USER_PROXY}} >> {0}".format(jdl_filename))
os.system("echo Arguments = {0} {1} {2} {3} {4} {5} {6} >> {7}".format( args.bar, args.firstRun, args.lastRun, args.biasVoltage, args.timeAlgo, tarball_name, outputDir, jdl_filename) )
os.system("echo Queue 1 >> {0}".format(jdl_filename))     
os.system("condor_submit {0}".format(jdl_filename))



# *** 4. Cleanup submission directory
print "\n##########     Cleanup submission directory     ##########\n"
os.system("rm *.jdl")
#os.system("mv {0} {1}/".format(tarball_name, args.outputDir) )
