# Estimation and histogram manipulation for multi-jet background


class MultijetEstimation(object):
    """Class for calculating multi-jet background"""
    def calcQCDsignal(self,histograms,name):
        """Calculate the QCD histogram"""
        ## 0.5*(yields["G"]/yields["A"] + yields["H"]/yields["B"])*yields["C"]
        h_A = histograms[name.format("0b0t")].Clone()
        h_B = histograms[name.format("1b0t")].Clone()
        h_C = histograms[name.format("2b0t")].Clone()
        h_G = histograms[name.format("0b2t")].Clone()
        h_H = histograms[name.format("1b2t")].Clone()

        qcd_hist = h_C.Clone()
        qcd_hist.Scale(0.5)

        h_G.Divide(h_A)
        h_H.Divide(h_B)

        h_G.Add(h_H)
        qcd_hist.Multiply(h_G)

        return qcd_hist

    def calcQCDvalidation(self,histograms,name):
        """Calculate the QCD histogram"""
        ## 0.5*(yields["D"]/yields["A"] + yields["E"]/yields["B"])*yields["C"]
        # print histograms.keys()
        h_A = histograms[name.format("0b0t")].Clone()
        h_B = histograms[name.format("1b0t")].Clone()
        h_C = histograms[name.format("2b0t")].Clone()
        h_D = histograms[name.format("0b1t")].Clone()
        h_E = histograms[name.format("1b1t")].Clone()

        qcd_hist = h_C.Clone()
        qcd_hist.Scale(0.5)

        h_D.Divide(h_A)
        h_E.Divide(h_B)

        h_D.Add(h_E)
        qcd_hist.Multiply(h_D)

        return qcd_hist
