"""Script to collect information about defects and CTI for individual ITL
sensors and use a metric to assign a score to sensors and rank them."""

from  eTraveler.clientAPI.connection import Connection
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import datetime
import collections
import numpy as np
from exploreRaft import exploreRaft

import argparse


class rank_ITL():

    def get_tests(self, htype = None, db = 'Prod', server = 'Prod', appSuffix = None):
        # Find all instances of the SR-EOT-02 traveler.

        if server == 'Prod': pS = True
        else: pS = False
        connect = Connection(operator='richard', db=db, exp='LSST-CAMERA', prodServer=pS)

        hardwareLabels = ['SR_Grade:', 'SR_Contract:']

        filePaths = connect.getFilepathsJH(htype='ITL-CCD', stepName='SR-CCD-EOT-05_Preflight-Check', travelerName='SR-CCD-EOT-05')#, hardwareLabels=hardwareLabels)

        # Defects
        # this step gives us bright columns and bright pixels,
        bright_defects_data  = connect.getResultsJH(htype='ITL-CCD', stepName = 'bright_defects_offline', travelerName='SR-EOT-02')
        # this step gives us dark columns and dark pixels
        dark_defects_data  = connect.getResultsJH(htype='ITL-CCD', stepName = 'dark_defects_offline', travelerName='SR-EOT-02')
        # this step gives us traps
        traps_data  = connect.getResultsJH(htype='ITL-CCD', stepName = 'traps_offline', travelerName='SR-EOT-02')

        # CTE
        # this step gives us cti_low_serial, cti_high_serial, cti_low_parallel, cti_high_parallel
        cte_data  = connect.getResultsJH(htype='ITL-CCD', stepName = 'cte_offline', travelerName='SR-EOT-02')

        # Read read_noise
        # this step gives us read_noise
        readnoise_data  = connect.getResultsJH(htype='ITL-CCD', stepName = 'read_noise_offline', travelerName='SR-EOT-02')

        nonlinearity_data = connect.getResultsJH(htype='ITL-CCD', stepName = 'flat_pairs_offline', travelerName='SR-EOT-02')

        ccd_list = []

        # Get a list of ccd's
        for ccd in filePaths:
            ccd_list.append(ccd)

        expDict = {}
        expDict['brightdefects']  = bright_defects_data
        expDict['darkdefects'] = dark_defects_data
        expDict['traps'] = traps_data
        expDict['cte'] = cte_data
        expDict['readnoise'] = readnoise_data
        expDict['nonlinearity'] = nonlinearity_data

        return ccd_list, expDict, filePaths

    def get_CTI(self, d, ccd):

        cti_low_serial = []
        cti_high_serial = []
        cti_low_parallel = []
        cti_high_parallel = []

        for amp in d['cte'][ccd]['steps']['cte_offline']['cte'][1:]:
            cti_low_serial.append(amp['cti_low_serial'])
            cti_high_serial.append(amp['cti_high_serial'])
            cti_low_parallel.append(amp['cti_low_parallel'])
            cti_high_parallel.append(amp['cti_high_parallel'])

        return np.array(cti_low_serial), np.array(cti_high_serial), np.array(cti_low_parallel), np.array(cti_high_parallel)

    def get_readnoise(self, d, ccd):

        read_noise = []

        for amp in d['readnoise'][ccd]['steps']['read_noise_offline']['read_noise'][1:]:
            read_noise.append(amp['read_noise'])

        return np.array(read_noise)

    def get_nonlinearity(self, d, ccd):

        nonlinearity = []

        for amp in d['nonlinearity'][ccd]['steps']['flat_pairs_offline']['flat_pairs'][1:]:
            nonlinearity.append(amp['max_frac_dev'])

        return np.array(nonlinearity)


    def get_bias(self, d, filePaths, ccd):
        # get the bias

        run = filePaths[ccd]['runInt']
        job_id = filePaths[ccd]['steps'][u'SR-CCD-EOT-05_Preflight-Check'][1]['activityId']
        vpath = filePaths[ccd]['steps'][u'SR-CCD-EOT-05_Preflight-Check'][12]['virtualPath'].strip('LSST/')
        summary_file = '/nfs/farm/g/lsst/u5/mirror/BNL-prod/prod/ITL-CCD/' + ccd + '/' +  str(run) + '/preflight_acq/v0/' + str(job_id) + '/summary.txt'
#        summary_file = '/Users/richard/lsst_archive/u5/mirror/BNL-prod/prod/ITL-CCD/' + ccd + '/' +  str(run) + '/preflight_acq/v0/' + str(job_id) + '/summary.txt'

        with open(summary_file) as f:
                lines = [line.split() for line in f]
        biases = np.array([float(line[6]) for line in lines[6:22]])

        bias = max(np.median(biases[0:8]),np.median(biases[8:]))

        return bias,biases


    def defectsFraction(self,d, ccd, col_len=2000, totPixels=1024000.0):
        #sum defect types and express as fraction (per segment)
        # totPixels = 1025024 E2V, 1024000 ITL

        amps = []
        numBrightPixels = []
        numBrightColumns = []
        numDarkPixels = []
        numDarkColumns = []
        numTraps = []

        for amp in d['brightdefects'][ccd]['steps']['bright_defects_offline']['bright_defects']:
            numBrightPixels.append(amp['bright_pixels'])
            numBrightColumns.append(amp['bright_columns'])

        for amp in d['darkdefects'][ccd]['steps']['dark_defects_offline']['dark_defects']:
            numDarkPixels.append(amp['dark_pixels'])
            numDarkColumns.append(amp['dark_columns'])

        for amp in d['traps'][ccd]['steps']['traps_offline']['traps']:
            numTraps.append(amp['num_traps'])

        numBrightPixels = np.array(numBrightPixels[1:])
        numBrightColumns = np.array(numBrightColumns[1:])
        numDarkPixels = np.array(numDarkPixels[1:])
        numDarkColumns = np.array(numDarkColumns[1:])
        numTraps = np.array(numTraps[1:])

        return ((numBrightPixels+numDarkPixels+numTraps+col_len*(numBrightColumns + numDarkColumns))/totPixels), numBrightColumns

    def metric(self,ccd_list, expDict, filePaths):

        eR = exploreRaft()

        #ccd_dic = {}
        RTM_dic={}
	bias_dic={}
        biases_dic={}
        for ccd in ccd_list:
            defectsfrac, nbcs = self.defectsFraction(expDict,ccd)
            nbc=nbcs.sum()

            ctils, ctihs, ctilp, ctihp = self.get_CTI(expDict,ccd)
            rn = self.get_readnoise(expDict,ccd)

            bias, biases = self.get_bias(expDict,filePaths, ccd)
            bias_dic[ccd]=bias
            biases_dic[ccd]=biases

            nonlinearity = self.get_nonlinearity(expDict, ccd)

            # Score = sCTI/2 + defects/0.2% + nonlinearity/1% + Bias/5000
            #Score = ctils.max()/2e-6 + defectsfrac.mean()/0.002 + max(nonlinearity)/.01 + bias/5000

            # Score = sCTI/2 + defectsfrac.mean()/0.002 + max(nonlineaerity)/.01% + bias/5000

            # Grading criteria: sensor is either SCIENCE (low-level noncompliance), RESERVE (compromised but usable in last resort), or ENGINEERING (do not use in FPA)
            # Based on typical ITL population where RN, sCTI, and defects are most common fails

            #if (rn.max()<10) and ((rn>8.5).sum()<3) and (1e6*ctils.max()<9) and (1e6*ctilp.max()<7) and (1e6*ctihs.max()<9) and (1e6*ctihp.max()<5) and (defectsfrac.mean()<0.049) and (nbc<5):
            #    GRADE="SCIENCE"
            #elif (rn.max()<14) and ((rn>11).sum()<4) and (1e6*ctils.max()<18) and ((1e6*ctils>10).sum()<4) and (1e6*ctilp.max()<10) and (1e6*ctihp.max()<10) and ((1e6*ctilp>10).sum()<4):
            #    GRADE="RESERVE"
            #else:
            #    GRADE="ENGIN."

            #ccd_dic[ccd] = GRADE
            parentRTM = eR.CCD_parent(ccd)
            #ccd_dic[ccd] = Score
            print ccd, parentRTM
            RTM_dic[ccd]= parentRTM

        #for key, value in sorted(ccd_dic.iteritems(), key=lambda (k,v): (v,k)):
            #print "%s: %s" % (key, value)

        return bias_dic, biases_dic, RTM_dic


if __name__ == "__main__":

    ## Command line arguments
    parser = argparse.ArgumentParser(
        description='Find archived data in the LSST  data Catalog. These include CCD test stand and vendor data files.')

    ##   The following are 'convenience options' which could also be specified in the filter string
    #parser.add_argument('-t', '--htype', default=None, help="hardware type (default=%(default)s)") #ITL-CCD
    parser.add_argument('-d', '--db', default='Prod', help="eT database (default=%(default)s)")
    parser.add_argument('-e', '--eTserver', default='Dev', help="eTraveler server (default=%(default)s)")
    parser.add_argument('--appSuffix', '--appSuffix', default='jrb',
                        help="eTraveler server (default=%(default)s)")
    args = parser.parse_args()

    rank = rank_ITL()
    ccd_list, expDict, filePaths = rank.get_tests(htype='ITL-CCD', db = args.db, server= args.eTserver, appSuffix = args.appSuffix)
    bias_dic, biases_dic, RTM_dic = rank.metric(ccd_list, expDict, filePaths)
    import pickle
    f=open('ITLnewdumpv2.p','wb')
    pickle.dump(expDict,f)
    pickle.dump(ccd_list,f)
    pickle.dump(bias_dic,f)
    pickle.dump(biases_dic,f)
    pickle.dump(RTM_dic,f)
f.close()
