#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 16:22:39 2021

@author: maria
"""
import numpy as np
import math

#Nup133
class SelectNup:
    
    def __init__(self, nup, term, model):
        
        
        refs = []
        ringmember = []
        self.model = model.upper()   
        
        if type(nup) == str: nup = tuple([nup])
        if type(term) == str: term = tuple([term])
        
        
        for i in range(len(nup)):
            tempnup = nup[i].upper()
            tempterm = term[i].upper()
            tempref = self.rotUnitCoords(tempnup, tempterm, self.model)
            ringmember.append(self.assignRing(tempref, tempnup, self.model))
            
            refs.append(tempref)    
        

        nrngs = [len(i) for i in refs] #number of rings per nup 
        refs = np.array(refs)
        
        self.nupindex = np.repeat(range(len(nrngs)), nrngs) #
        
        self.refs = np.concatenate(refs)
        self.ringmember = np.concatenate(ringmember)
        

        #self.refs = np.reshape(refs, [rshape[0]*rshape[1], rshape[2]])
        self.z = self.distz(self.refs)
        self.r = self.radii(self.refs)

        self.ringAngles = self.allang(self.refs)
        self.rotAng, self.rotOffset = self.rotang(self.z, self.ringAngles)
     
   
    def rotUnitCoords(self, nup, term, model):
        #MODEL 5A9Q
        #Nup43 N getcr #1.1/0,R,9,I:4@CA
        
        if model == "SMILE":
            if nup == "NOTANUP":
                if term == "NOTATERM":
                    lefteye = np.array([3250, 2500, 0]) 
                    righteye = np.array([4750, 2500, 0]) 
                    nose1 = np.array([4000, 1750, 0]) 
                    tooth1 = np.array([3000, 1500, 0]) 
                    tooth2 = np.array([3500, 1150, 0]) 
                    tooth3 = np.array([4000, 1000, 0]) 
                    tooth4 = np.array([4500, 1150, 0]) 
                    tooth5 = np.array([5000, 1500, 0]) 
                    ref = np.array([lefteye, righteye, nose1, tooth1, tooth2, tooth3, tooth4, tooth5])
                    
        if (model == "SIMPLE"): # Simple model for implementation purposes
            if (type(nup) == str):
                if (type(term) == str):
                    auth = np.array([1196.307, 578.280, 1])
                                        
                    ref = np.array([auth])

        if (model == "OTHER"): # Simple model for implementation purposes
            if (type(nup) == tuple):
                if (type(term) == tuple):
                    auth1 = np.array([1., 500., 1.])
                    auth2 = np.array([50., 600., 1.])
                    auth3 = np.array([80., 450., 1.])
                    auth4 = np.array([140., 550., 1.])
                    auth5 = np.array([190., 400., 1.])
                    auth6 = np.array([185., 950., 1.])
                    ref = np.array([auth1, auth2, auth3, auth4, auth5, auth6])
            
            
            
        if(model == "7R5K"): #alpha-fold, human NPC, constricted
            if (nup == "RANBP2"):
            	if (term == "N"):
            		authN_00 = np.array([510.043,832.733,1231.634])
            		authN_01 = np.array([445.808,925.765,1258.438])
            		authN_02 = np.array([569.931,840.576,1264.304])
            		authN_03 = np.array([503.951,941.958,1274.828])
            		authN_04 = np.array([428.437,849.196,1293.412])
            		ref = np.array([authN_00, authN_01, authN_02, authN_03, authN_04])
            
            	if (term == "C"):
            		authC_00 = np.array([504.813,743.724,1303.004])
            		authC_01 = np.array([424.021,836.440,1310.132])
            		authC_02 = np.array([491.417,798.700,1329.034])
            		authC_03 = np.array([459.114,872.021,1343.155])
            		authC_04 = np.array([527.402,838.200,1344.866])
            		ref = np.array([authC_00, authC_01, authC_02, authC_03, authC_04])
            
            if (nup == "NUP210"):
            	if (term == "N"):
            		authN_10 = np.array([370.788,1262.682,968.394])
            		authN_11 = np.array([343.940,1209.459,959.211])
            		authN_12 = np.array([329.215,1152.709,952.039])
            		authN_13 = np.array([317.326,1092.081,938.016])
            		authN_14 = np.array([317.540,838.407,1012.491])
            		authN_15 = np.array([331.430,778.375,1007.137])
            		authN_16 = np.array([348.014,721.824,999.505])
            		authN_17 = np.array([378.970,670.048,984.231])
            		ref = np.array([authN_10, authN_11, authN_12, authN_13, authN_14, authN_15, authN_16, authN_17])
            
            	if (term == "C"):
            		authC_10 = np.array([578.211,1057.786,980.116])
            		authC_11 = np.array([565.344,1030.077,972.762])
            		authC_12 = np.array([554.804,996.797,993.606])
            		authC_13 = np.array([550.638,982.269,1005.711])
            		authC_14 = np.array([555.152,958.843,959.020])
            		authC_15 = np.array([554.762,940.445,946.648])
            		authC_16 = np.array([570.751,910.371,984.501])
            		authC_17 = np.array([585.630,884.535,987.569])
            		ref = np.array([authC_10, authC_11, authC_12, authC_13, authC_14, authC_15, authC_16, authC_17])
            
            if (nup == "ALADIN"):
            	if (term == "N"):
            		authN_40 = np.array([590.792,1068.107,908.350])
            		authN_41 = np.array([594.803,864.042,1052.684])
            		ref = np.array([authN_40, authN_41])
            
            	if (term == "C"):
            		authC_40 = np.array([599.106,1049.022,911.680])
            		authC_41 = np.array([603.062,884.550,1049.119])
            		ref = np.array([authC_40, authC_41])
            
            
            if (nup == "NUP93"):
            	if (term == "N"):
                    authN_A0 = np.array([693.120,998.544,998.190])
                    authN_A1 = np.array([705.093,926.668,1027.205])
                    authN_A2 = np.array([691.636,938.949,965.499])
                    authN_A3 = np.array([704.099,1015.256,931.436])
                    authN_A4 = np.array([570.953,985.442,1193.802])
                    authN_A5 = np.array([544.261,1064.634,1243.666])
                    authN_A6 = np.array([549.040,853.400,735.036])
                    ref = np.array([authN_A0, authN_A1, authN_A2, authN_A3, authN_A4, authN_A5, authN_A6])
            
            	if (term == "C"):
                    authC_A0 = np.array([622.047,1018.485,994.302])
                    authC_A1 = np.array([619.321,928.973,1064.974])
                    authC_A2 = np.array([618.512,926.555,976.617])
                    authC_A3 = np.array([614.854,1011.768,893.976])
                    authC_A4 = np.array([495.450,1119.245,1212.592])
                    authC_A5 = np.array([495.681,1210.853,1240.324])
                    authC_A6 = np.array([515.746,713.244,741.804])
                    ref = np.array([authC_A0, authC_A1, authC_A2, authC_A3, authC_A4, authC_A5, authC_A6])

            if (nup == "NUP188"):
            	if (term == "N"):
            		authN_B0 = np.array([709.756,980.766,1069.307])
            		authN_B1 = np.array([705.844,958.936,892.389])
            		ref = np.array([authN_B0, authN_B1])
            
            	if (term == "C"):
            		authC_B0 = np.array([681.105,857.524,1056.848])
            		authC_B1 = np.array([676.423,1082.075,900.732])
            		ref = np.array([authC_B0, authC_B1])
            
            if (nup == "NUP205"):
            	if (term == "N"):
            		authN_C0 = np.array([710.152,1008.472,1046.729])
            		authN_C1 = np.array([708.160,928.333,915.040])
            		authN_C2 = np.array([538.763,962.018,1155.895])
            		authN_C3 = np.array([523.299,1065.300,1193.237])
            		authN_C4 = np.array([523.804,852.209,783.399])
            		ref = np.array([authN_C0, authN_C1, authN_C2, authN_C3, authN_C4])
            
            	if (term == "C"):
            		authC_C0 = np.array([633.188,928.264,993.209])
            		authC_C1 = np.array([633.415,1009.031,970.933])
            		authC_C2 = np.array([575.332,1074.750,1199.786])
            		authC_C3 = np.array([586.635,1137.907,1274.827])
            		authC_C4 = np.array([597.361,782.751,707.561])
            		ref = np.array([authC_C0, authC_C1, authC_C2, authC_C3, authC_C4])
            
            if (nup == "NUP155"):
            	if (term == "N"):
            		authN_D0 = np.array([647.172,1062.852,973.846])
            		authN_D1 = np.array([588.111,968.266,1053.434])
            		authN_D2 = np.array([643.592,874.960,981.972])
            		authN_D3 = np.array([586.063,959.949,909.277])
            		authN_D4 = np.array([559.008,966.031,1125.355])
            		authN_D5 = np.array([574.234,979.147,849.650])
            		ref = np.array([authN_D0, authN_D1, authN_D2, authN_D3, authN_D4, authN_D5])
            
            	if (term == "C"):
            		authC_D0 = np.array([627.762,975.219,1043.680])
            		authC_D1 = np.array([663.624,902.004,1097.449])
            		authC_D2 = np.array([619.624,967.270,919.613])
            		authC_D3 = np.array([656.588,1038.417,861.813])
            		authC_D4 = np.array([501.048,931.450,1202.022])
            		authC_D5 = np.array([491.287,987.465,766.724])
            		ref = np.array([authC_D0, authC_D1, authC_D2, authC_D3, authC_D4, authC_D5])
            
            if (nup == "NDL1"):
            	if (term == "N"):
            		authN_E0 = np.array([563.694,1044.522,985.531])
            		authN_E1 = np.array([551.410,896.729,986.573])
            		ref = np.array([authN_E0, authN_E1])
            
            	if (term == "C"):
            		authC_E0 = np.array([578.011,1028.859,924.543])
            		authC_E1 = np.array([580.626,906.410,1043.326])
            		ref = np.array([authC_E0, authC_E1])
            
            if (nup == "NUP35"):
            	if (term == "N"):
            		authN_F0 = np.array([676.528,1105.014,899.746])
            		authN_F1 = np.array([637.928,1028.563,932.945])
            		authN_F2 = np.array([682.469,835.707,1063.081])
            		authN_F3 = np.array([637.858,911.216,1024.657])
            		ref = np.array([authN_F0, authN_F1, authN_F2, authN_F3])
            
            	if (term == "C"):
            		authC_F0 = np.array([559.777,1071.516,920.964])
            		authC_F1 = np.array([565.726,996.840,969.631])
            		authC_F2 = np.array([561.596,863.885,1047.702])
            		authC_F3 = np.array([560.477,952.193,997.591])
            		ref = np.array([authC_F0, authC_F1, authC_F2, authC_F3])
            
            if (nup == "NUP54"):
            	if (term == "N"):
            		authN_H0 = np.array([751.321,979.227,985.578])
            		authN_H1 = np.array([742.708,894.937,988.245])
            		authN_H2 = np.array([749.783,959.981,974.928])
            		authN_H3 = np.array([738.824,1048.293,972.078])
            		ref = np.array([authN_H0, authN_H1, authN_H2, authN_H3])
            
            	if (term == "C"):
            		authC_H0 = np.array([692.384,1027.067,952.831])
            		authC_H1 = np.array([708.593,970.004,996.397])
            		authC_H2 = np.array([693.771,910.933,1011.087])
            		authC_H3 = np.array([707.296,972.150,962.872])
            		ref = np.array([authC_H0, authC_H1, authC_H2, authC_H3])
            
            if (nup == "NUP58"):
            	if (term == "N"):
            		authN_I0 = np.array([748.190,968.367,1018.800])
            		authN_I1 = np.array([741.904,872.797,1015.517])
            		authN_I2 = np.array([744.592,970.322,941.741])
            		authN_I3 = np.array([738.783,1070.108,944.618])
            		ref = np.array([authN_I0, authN_I1, authN_I2, authN_I3])
            
            	if (term == "C"):
            		authC_I0 = np.array([696.619,1016.209,949.657])
            		authC_I1 = np.array([706.432,960.917,988.972])
            		authC_I2 = np.array([697.906,921.848,1013.839])
            		authC_I3 = np.array([704.406,981.262,970.180])
            		ref = np.array([authC_I0, authC_I1, authC_I2, authC_I3])
            
            if (nup == "NUP62"):
            	if (term == "N"):
            		authN_J0 = np.array([737.974,969.303,1020.178])
            		authN_J1 = np.array([733.801,876.635,1020.709])
            		authN_J2 = np.array([734.309,969.116,941.007])
            		authN_J3 = np.array([731.112,1066.058,939.030])
            		authN_J4 = np.array([636.202,998.038,1276.585])
            		ref = np.array([authN_J0, authN_J1, authN_J2, authN_J3, authN_J4])
            
            	if (term == "C"):
            		authC_J0 = np.array([683.679,1019.575,950.773])
            		authC_J1 = np.array([697.269,968.076,995.719])
            		authC_J2 = np.array([685.008,918.067,1013.550])
            		authC_J3 = np.array([695.837,973.867,962.994])
            		authC_J4 = np.array([584.290,1064.969,1260.496])
            		ref = np.array([authC_J0, authC_J1, authC_J2, authC_J3, authC_J4])
            
            if (nup == "NUP133"):
            	if (term == "N"):
            		authN_K0 = np.array([569.935,615.923,1162.184])
            		authN_K1 = np.array([482.969,669.808,1178.008])
            		authN_K2 = np.array([571.414,1282.931,811.536])
            		authN_K3 = np.array([483.200,1245.573,792.518])
            		ref = np.array([authN_K0, authN_K1, authN_K2, authN_K3])
            
            	if (term == "C"):
            		authC_K0 = np.array([479.978,727.684,1209.271])
            		authC_K1 = np.array([407.932,837.562,1242.624])
            		authC_K2 = np.array([458.660,1187.939,757.830])
            		authC_K3 = np.array([403.947,1075.697,732.736])
            		ref = np.array([authC_K0, authC_K1, authC_K2, authC_K3])
            
            if (nup == "NUP107"):
            	if (term == "N"):
            		authN_L0 = np.array([559.908,768.369,1246.970])
            		authN_L1 = np.array([494.829,866.444,1279.050])
            		authN_L2 = np.array([545.927,1149.778,731.379])
            		authN_L3 = np.array([487.595,1051.594,697.217])
            		ref = np.array([authN_L0, authN_L1, authN_L2, authN_L3])
            
            	if (term == "C"):
            		authC_L0 = np.array([519.183,680.358,1222.889])
            		authC_L1 = np.array([438.427,783.495,1254.615])
            		authC_L2 = np.array([500.564,1233.199,751.775])
            		authC_L3 = np.array([430.825,1132.817,720.023])
            		ref = np.array([authC_L0, authC_L1, authC_L2, authC_L3])
            
            if (nup == "NUP96"):
            	if (term == "N"):
            		authN_M0 = np.array([552.847,914.404,1245.745])
            		authN_M1 = np.array([494.824,1010.898,1272.903])
            		authN_M2 = np.array([552.038,1000.567,736.947])
            		authN_M3 = np.array([497.981,905.351,702.802])
            		ref = np.array([authN_M0, authN_M1, authN_M2, authN_M3])
            
            	if (term == "C"):
            		authC_M0 = np.array([452.004,932.756,1215.370])
            		authC_M1 = np.array([405.398,1031.436,1235.971])
            		authC_M2 = np.array([453.649,989.636,765.209])
            		authC_M3 = np.array([407.457,883.188,734.412])
            		ref = np.array([authC_M0, authC_M1, authC_M2, authC_M3])
            
            if (nup == "SEC13"):
            	if (term == "N"):
            		authN_N0 = np.array([537.184,882.463,1247.806])
            		authN_N1 = np.array([477.071,980.306,1274.913])
            		authN_N2 = np.array([536.969,1032.591,733.808])
            		authN_N3 = np.array([478.872,935.058,699.590])
            		ref = np.array([authN_N0, authN_N1, authN_N2, authN_N3])
            
            	if (term == "C"):
            		authC_N0 = np.array([543.846,879.698,1256.947])
            		authC_N1 = np.array([483.231,976.960,1284.051])
            		authC_N2 = np.array([544.676,1035.317,725.750])
            		authC_N3 = np.array([486.030,938.986,691.331])
            		ref = np.array([authC_N0, authC_N1, authC_N2, authC_N3])
            
            if (nup == "SEH1"):
            	if (term == "N"):
            		authN_O0 = np.array([557.582,977.776,1282.196])
            		authN_O1 = np.array([499.989,1075.037,1312.512])
            		authN_O2 = np.array([562.810,940.266,702.330])
            		authN_O3 = np.array([512.269,844.769,665.310])
            		ref = np.array([authN_O0, authN_O1, authN_O2, authN_O3])
            
            	if (term == "C"):
            		authC_O0 = np.array([555.403,973.760,1298.380])
            		authC_O1 = np.array([495.400,1072.448,1329.595])
            		authC_O2 = np.array([562.840,945.656,684.978])
            		authC_O3 = np.array([510.680,849.362,648.079])
            		ref = np.array([authC_O0, authC_O1, authC_O2, authC_O3])
            
            if (nup == "NUP85"):
            	if (term == "N"):
            		authN_P0 = np.array([580.226,974.519,1271.114])
            		authN_P1 = np.array([522.213,1067.511,1305.033])
            		authN_P2 = np.array([583.164,944.567,715.615])
            		authN_P3 = np.array([532.896,851.691,677.210])
            		ref = np.array([authN_P0, authN_P1, authN_P2, authN_P3])
            
            	if (term == "C"):
            		authC_P0 = np.array([495.807,954.611,1210.612])
            		authC_P1 = np.array([442.613,1053.953,1237.425])
            		authC_P2 = np.array([493.911,960.633,764.306])
            		authC_P3 = np.array([445.699,858.887,731.979])
            		ref = np.array([authC_P0, authC_P1, authC_P2, authC_P3])
            
            if (nup == "NUP43"):
            	if (term == "N"):
            		authN_Q0 = np.array([545.254,1014.513,1255.226])
            		authN_Q1 = np.array([495.656,1109.654,1282.469])
            		authN_Q2 = np.array([547.928,902.715,725.079])
            		authN_Q3 = np.array([504.418,806.473,689.787])
            		ref = np.array([authN_Q0, authN_Q1, authN_Q2, authN_Q3])
            
            	if (term == "C"):
            		authC_Q0 = np.array([539.937,1017.560,1252.388])
            		authC_Q1 = np.array([490.688,1113.342,1278.690])
            		authC_Q2 = np.array([543.404,899.819,726.995])
            		authC_Q3 = np.array([499.279,802.304,692.154])
            		ref = np.array([authC_Q0, authC_Q1, authC_Q2, authC_Q3])
            
            if (nup == "NUP160"):
            	if (term == "N"):
            		authN_R0 = np.array([433.129,946.632,1159.874])
            		authN_R1 = np.array([393.069,1048.780,1177.163])
            		authN_R2 = np.array([432.113,962.600,816.734])
            		authN_R3 = np.array([393.422,861.689,786.505])
            		ref = np.array([authN_R0, authN_R1, authN_R2, authN_R3])
            
            	if (term == "C"):
            		authC_R0 = np.array([532.120,905.516,1290.860])
            		authC_R1 = np.array([473.016,1003.872,1317.762])
            		authC_R2 = np.array([536.174,1010.546,690.125])
            		authC_R3 = np.array([481.344,912.989,655.786])
            		ref = np.array([authC_R0, authC_R1, authC_R2, authC_R3])
            
            if (nup == "NUP37"):
            	if (term == "N"):
            		authN_S0 = np.array([491.973,965.842,1188.629])
            		authN_S1 = np.array([443.888,1067.185,1219.103])
            		authN_S2 = np.array([491.543,951.063,777.149])
            		authN_S3 = np.array([448.624,848.874,748.916])
            		ref = np.array([authN_S0, authN_S1, authN_S2, authN_S3])
            
            	if (term == "C"):
            		authC_S0 = np.array([487.411,984.669,1191.577])
            		authC_S1 = np.array([440.538,1086.611,1218.938])
            		authC_S2 = np.array([490.137,931.696,776.153])
            		authC_S3 = np.array([448.242,829.200,750.709])
            		ref = np.array([authC_S0, authC_S1, authC_S2, authC_S3])
            
            if (nup == "ELYS"):
            	if (term == "N"):
            		authN_T0 = np.array([405.053,910.997,811.798])
            		authN_T1 = np.array([368.835,808.965,798.442])
            		ref = np.array([authN_T0, authN_T1])
            
            	if (term == "C"):
            		authC_T0 = np.array([479.992,979.693,724.077])
            		authC_T1 = np.array([428.933,868.979,693.913])
            		ref = np.array([authC_T0, authC_T1])
            
            if (nup == "NUP98"):
            	if (term == "N"):
            		authN_U0 = np.array([658.239,1001.835,1198.455])
            		authN_U1 = np.array([500.401,956.959,1158.438])
            		authN_U2 = np.array([642.058,1022.906,1021.722])
            		authN_U3 = np.array([636.836,945.966,1079.213])
            		authN_U4 = np.array([636.411,918.901,938.572])
            		authN_U5 = np.array([630.639,994.096,880.254])
            		authN_U6 = np.array([516.367,965.292,809.463])
            		ref = np.array([authN_U0, authN_U1, authN_U2, authN_U3, authN_U4, authN_U5, authN_U6])
            
            	if (term == "C"):
            		authC_U0 = np.array([649.104,999.358,1223.761])
            		authC_U1 = np.array([489.426,957.062,1131.646])
            		authC_U2 = np.array([628.084,1046.210,1031.725])
            		authC_U3 = np.array([619.898,959.424,1098.433])
            		authC_U4 = np.array([624.459,895.561,926.343])
            		authC_U5 = np.array([614.524,979.810,860.898])
            		authC_U6 = np.array([501.532,965.166,834.328])
            		ref = np.array([authC_U0, authC_U1, authC_U2, authC_U3, authC_U4, authC_U5, authC_U6])
            
            if (nup == "NUP214"):
            	if (term == "N"):
            		authN_V0 = np.array([653.506,998.119,1281.827])
            		ref = np.array([authN_V0])
            
            	if (term == "C"):
            		authC_V0 = np.array([649.416,1011.751,1281.999])
            		ref = np.array([authC_V0])
            
            if (nup == "NUP88"):
            	if (term == "N"):
            		authN_W0 = np.array([604.570,1012.598,1276.027])
            		ref = np.array([authN_W0])
            
            	if (term == "C"):
            		authC_W0 = np.array([589.418,1064.393,1269.213])
            		ref = np.array([authC_W0])
                
            
        if(model == "7R5J"):    # alpha-fold model, human NPC, dilated
            if (nup == "RANBP2"):
            	if (term == "N"):
            		authN_00 = np.array([484.356,991.167,1249.945])
            		authN_01 = np.array([446.498,1097.835,1283.338])
            		authN_02 = np.array([549.894,990.023,1276.608])
            		authN_03 = np.array([513.609,1104.964,1286.513])
            		authN_04 = np.array([418.023,1030.990,1316.641])
            		ref = np.array([authN_00, authN_01, authN_02, authN_03, authN_04])
            
            	if (term == "C"):
            		authC_00 = np.array([467.568,911.839,1321.510])
            		authC_01 = np.array([409.783,1020.856,1334.846])
            		authC_02 = np.array([476.682,965.246,1364.860])
            		authC_03 = np.array([456.609,1047.712,1363.927])
            		authC_04 = np.array([514.315,992.783,1362.120])
            		ref = np.array([authC_00, authC_01, authC_02, authC_03, authC_04])
            
            if (nup == "NUP210"):
            	if (term == "N"):
            		authN_10 = np.array([483.138,1512.950,952.884])
            		authN_11 = np.array([436.905,1465.193,942.890])
            		authN_12 = np.array([408.774,1409.150,936.054])
            		authN_13 = np.array([397.516,1352.999,934.545])
            		authN_14 = np.array([298.193,1016.994,1072.259])
            		authN_15 = np.array([289.578,953.466,1064.153])
            		authN_16 = np.array([279.562,899.424,1058.211])
            		authN_17 = np.array([265.462,833.896,1047.998])
            		ref = np.array([authN_10, authN_11, authN_12, authN_13, authN_14, authN_15, authN_16, authN_17])
            
            	if (term == "C"):
            		authC_10 = np.array([575.245,1203.292,947.364])
            		authC_11 = np.array([553.706,1186.357,982.139])
            		authC_12 = np.array([534.240,1154.949,1000.082])
            		authC_13 = np.array([517.777,1131.089,1014.099])
            		authC_14 = np.array([516.332,1136.061,926.119])
            		authC_15 = np.array([518.193,1111.818,932.473])
            		authC_16 = np.array([516.786,1074.505,953.731])
            		authC_17 = np.array([532.806,1051.182,997.704])
            		ref = np.array([authC_10, authC_11, authC_12, authC_13, authC_14, authC_15, authC_16, authC_17])
            
            if (nup == "ALADIN"):
            	if (term == "N"):
            		authN_40 = np.array([597.520,1210.312,901.557])
            		authN_41 = np.array([540.978,1003.410,1049.895])
            		ref = np.array([authN_40, authN_41])
            
            	if (term == "C"):
            		authC_40 = np.array([602.008,1189.651,905.846])
            		authC_41 = np.array([555.524,1019.406,1047.318])
            		ref = np.array([authC_40, authC_41])
            
            if (nup == "NUP93"):
            	if (term == "N"):
            		authN_A0 = np.array([674.420,1111.901,990.369])
            		authN_A1 = np.array([668.557,1034.423,1023.867])
            		authN_A2 = np.array([664.295,1053.451,967.036])
            		authN_A3 = np.array([692.481,1127.041,932.147])
            		authN_A4 = np.array([580.057,1126.571,1200.047])
            		authN_A5 = np.array([581.430,1210.147,1250.213])
            		authN_A6 = np.array([531.399,973.428,721.050])
            		ref = np.array([authN_A0, authN_A1, authN_A2, authN_A3, authN_A4, authN_A5, authN_A6])
            
            	if (term == "C"):
            		authC_A0 = np.array([609.711,1147.915,994.735])
            		authC_A1 = np.array([580.379,1058.881,1065.389])
            		authC_A2 = np.array([585.551,1058.953,967.315])
            		authC_A3 = np.array([607.812,1144.050,886.914])
            		authC_A4 = np.array([553.198,1270.773,1231.144])
            		authC_A5 = np.array([591.669,1377.244,1265.375])
            		authC_A6 = np.array([458.816,825.014,696.912])
            		ref = np.array([authC_A0, authC_A1, authC_A2, authC_A3, authC_A4, authC_A5, authC_A6])
            
            if (nup == "NUP188"):
            	if (term == "N"):
            		authN_B0 = np.array([683.475,1086.125,1067.435])
            		authN_B1 = np.array([682.765,1072.513,889.923])
            		ref = np.array([authN_B0, authN_B1])
            
            	if (term == "C"):
            		authC_B0 = np.array([626.389,974.613,1048.497])
            		authC_B1 = np.array([683.878,1198.532,902.348])
            		ref = np.array([authC_B0, authC_B1])
            
            if (nup == "NUP205"):
            	if (term == "N"):
            		authN_C0 = np.array([691.739,1114.212,1046.153])
            		authN_C1 = np.array([676.158,1041.000,911.395])
            		authN_C2 = np.array([538.643,1112.686,1166.726])
            		authN_C3 = np.array([555.344,1216.054,1202.981])
            		authN_C4 = np.array([502.441,979.661,766.892])
            		ref = np.array([authN_C0, authN_C1, authN_C2, authN_C3, authN_C4])
            
            	if (term == "C"):
            		authC_C0 = np.array([600.481,1056.835,985.969])
            		authC_C1 = np.array([620.956,1136.054,967.444])
            		authC_C2 = np.array([608.645,1211.414,1203.102])
            		authC_C3 = np.array([644.917,1268.926,1274.860])
            		authC_C4 = np.array([567.192,893.298,701.025])
            		ref = np.array([authC_C0, authC_C1, authC_C2, authC_C3, authC_C4])
            
            if (nup == "NUP155"):
            	if (term == "N"):
            		authN_D0 = np.array([650.562,1185.234,971.280])
            		authN_D1 = np.array([565.453,1112.784,1052.431])
            		authN_D2 = np.array([599.907,999.455,975.478])
            		authN_D3 = np.array([554.257,1099.887,900.577])
            		authN_D4 = np.array([542.445,1107.597,1123.531])
            		authN_D5 = np.array([546.740,1101.349,847.779])
            		ref = np.array([authN_D0, authN_D1, authN_D2, authN_D3, authN_D4, authN_D5])
            
            	if (term == "C"):
            		authC_D0 = np.array([594.628,1107.277,1038.985])
            		authC_D1 = np.array([612.347,1016.131,1087.060])
            		authC_D2 = np.array([596.271,1099.646,914.790])
            		authC_D3 = np.array([653.302,1165.453,852.727])
            		authC_D4 = np.array([495.034,1085.503,1219.276])
            		authC_D5 = np.array([485.946,1118.453,746.177])
            		ref = np.array([authC_D0, authC_D1, authC_D2, authC_D3, authC_D4, authC_D5])
            
            if (nup == "NDL1"):
            	if (term == "N"):
            		authN_E0 = np.array([559.647,1193.138,976.618])
            		authN_E1 = np.array([519.633,1042.564,976.723])
            		ref = np.array([authN_E0, authN_E1])
            
            	if (term == "C"):
            		authC_E0 = np.array([576.756,1173.104,917.608])
            		authC_E1 = np.array([542.252,1047.085,1037.103])
            		ref = np.array([authC_E0, authC_E1])
            
            if (nup == "NUP35"):
            	if (term == "N"):
            		authN_F0 = np.array([683.715,1225.869,902.303])
            		authN_F1 = np.array([630.132,1154.104,926.594])
            		authN_F2 = np.array([623.577,955.359,1054.934])
            		authN_F3 = np.array([597.796,1037.611,1014.954])
            		ref = np.array([authN_F0, authN_F1, authN_F2, authN_F3])
            
            	if (term == "C"):
            		authC_F0 = np.array([566.189,1219.138,911.362])
            		authC_F1 = np.array([557.522,1132.797,963.747])
            		authC_F2 = np.array([510.034,1012.200,1039.707])
            		authC_F3 = np.array([535.276,1092.908,990.309])
            		ref = np.array([authC_F0, authC_F1, authC_F2, authC_F3])
            
            if (nup == "NUP54"):
            	if (term == "N"):
            		authN_H0 = np.array([724.867,1075.495,981.602])
            		authN_H1 = np.array([696.967,997.572,981.446])
            		authN_H2 = np.array([726.140,1053.081,978.408])
            		authN_H3 = np.array([732.052,1148.241,976.183])
            		ref = np.array([authN_H0, authN_H1, authN_H2, authN_H3])
            
            	if (term == "C"):
            		authC_H0 = np.array([690.816,1144.626,951.525])
            		authC_H1 = np.array([680.179,1077.874,995.470])
            		authC_H2 = np.array([656.712,1016.532,1005.113])
            		authC_H3 = np.array([682.135,1083.859,961.439])
            		ref = np.array([authC_H0, authC_H1, authC_H2, authC_H3])
            
            if (nup == "NUP58"):
            	if (term == "N"):
            		authN_I0 = np.array([712.425,1062.305,1012.030])
            		authN_I1 = np.array([692.492,974.928,1007.658])
            		authN_I2 = np.array([725.150,1071.570,948.285])
            		authN_I3 = np.array([739.605,1170.164,949.493])
            		ref = np.array([authN_I0, authN_I1, authN_I2, authN_I3])
            
            	if (term == "C"):
            		authC_I0 = np.array([692.236,1133.896,947.666])
            		authC_I1 = np.array([676.385,1070.000,987.902])
            		authC_I2 = np.array([664.060,1024.427,1009.336])
            		authC_I3 = np.array([681.875,1092.853,968.687])
            		ref = np.array([authC_I0, authC_I1, authC_I2, authC_I3])
            
            if (nup == "NUP62"):
            	if (term == "N"):
            		authN_J0 = np.array([702.912,1066.712,1011.850])
            		authN_J1 = np.array([685.613,980.297,1013.304])
            		authN_J2 = np.array([714.888,1073.601,948.011])
            		authN_J3 = np.array([731.324,1168.207,943.484])
            		authN_J4 = np.array([647.165,1125.538,1257.712])
            		ref = np.array([authN_J0, authN_J1, authN_J2, authN_J3, authN_J4])
            
            	if (term == "C"):
            		authC_J0 = np.array([680.983,1142.024,947.225])
            		authC_J1 = np.array([669.044,1079.250,995.513])
            		authC_J2 = np.array([650.130,1024.431,1009.233])
            		authC_J3 = np.array([671.680,1087.560,960.728])
            		authC_J4 = np.array([609.460,1203.067,1260.504])
            		ref = np.array([authC_J0, authC_J1, authC_J2, authC_J3, authC_J4])
            
            if (nup == "NUP133"):
            	if (term == "N"):
            		authN_K0 = np.array([485.589,750.478,1188.752])
            		authN_K1 = np.array([395.898,827.652,1202.779])
            		authN_K2 = np.array([646.994,1419.842,803.587])
            		authN_K3 = np.array([536.572,1389.385,764.993])
            		ref = np.array([authN_K0, authN_K1, authN_K2, authN_K3])
            
            	if (term == "C"):
            		authC_K0 = np.array([425.847,902.375,1231.626])
            		authC_K1 = np.array([390.521,1026.607,1270.424])
            		authC_K2 = np.array([508.220,1322.842,731.172])
            		authC_K3 = np.array([430.798,1220.084,692.144])
            		ref = np.array([authC_K0, authC_K1, authC_K2, authC_K3])
            
            if (nup == "NUP107"):
            	if (term == "N"):
            		authN_L0 = np.array([520.991,919.384,1252.135])
            		authN_L1 = np.array([487.138,1031.566,1297.824])
            		authN_L2 = np.array([590.643,1256.501,701.266])
            		authN_L3 = np.array([513.243,1176.557,670.501])
            		ref = np.array([authN_L0, authN_L1, authN_L2, authN_L3])
            
            	if (term == "C"):
            		authC_L0 = np.array([451.640,846.821,1243.413])
            		authC_L1 = np.array([407.134,966.744,1280.932])
            		authC_L2 = np.array([559.721,1350.474,724.695])
            		authC_L3 = np.array([471.791,1269.143,684.881])
            		ref = np.array([authC_L0, authC_L1, authC_L2, authC_L3])
            
            if (nup == "NUP96"):
            	if (term == "N"):
            		authN_M0 = np.array([550.097,1061.035,1256.776])
            		authN_M1 = np.array([521.570,1169.746,1289.101])
            		authN_M2 = np.array([555.811,1114.817,730.993])
            		authN_M3 = np.array([494.446,1033.731,680.401])
            		ref = np.array([authN_M0, authN_M1, authN_M2, authN_M3])
            
            	if (term == "C"):
            		authC_M0 = np.array([461.291,1104.500,1237.538])
            		authC_M1 = np.array([438.050,1222.926,1269.250])
            		authC_M2 = np.array([455.175,1112.518,721.235])
            		authC_M3 = np.array([396.180,1015.676,694.460])
            		ref = np.array([authC_M0, authC_M1, authC_M2, authC_M3])
            
            if (nup == "SEC13"):
            	if (term == "N"):
            		authN_N0 = np.array([528.012,1034.587,1260.349])
            		authN_N1 = np.array([496.676,1146.215,1294.160])
            		authN_N2 = np.array([548.970,1146.785,719.478])
            		authN_N3 = np.array([480.576,1065.093,675.430])
            		ref = np.array([authN_N0, authN_N1, authN_N2, authN_N3])
            
            	if (term == "C"):
            		authC_N0 = np.array([533.439,1030.252,1269.004])
            		authC_N1 = np.array([501.690,1141.495,1302.874])
            		authC_N2 = np.array([557.968,1147.351,713.058])
            		authC_N3 = np.array([487.493,1067.842,667.179])
            		ref = np.array([authC_N0, authC_N1, authC_N2, authC_N3])
            
            if (nup == "SEH1"):
            	if (term == "N"):
            		authN_O0 = np.array([575.148,1118.599,1292.904])
            		authN_O1 = np.array([549.148,1232.768,1325.869])
            		authN_O2 = np.array([563.497,1056.386,687.057])
            		authN_O3 = np.array([501.001,966.361,644.267])
            		ref = np.array([authN_O0, authN_O1, authN_O2, authN_O3])
            
            	if (term == "C"):
            		authC_O0 = np.array([574.051,1112.151,1308.617])
            		authC_O1 = np.array([545.290,1229.808,1343.008])
            		authC_O2 = np.array([565.801,1060.884,669.624])
            		authC_O3 = np.array([502.018,971.779,627.387])
            		ref = np.array([authC_O0, authC_O1, authC_O2, authC_O3])
            
            if (nup == "NUP85"):
            	if (term == "N"):
            		authN_P0 = np.array([594.386,1113.270,1278.046])
            		authN_P1 = np.array([568.427,1221.085,1315.710])
            		authN_P2 = np.array([583.040,1058.507,702.000])
            		authN_P3 = np.array([521.593,968.352,658.497])
            		ref = np.array([authN_P0, authN_P1, authN_P2, authN_P3])
            
            	if (term == "C"):
            		authC_P0 = np.array([505.623,1114.455,1226.057])
            		authC_P1 = np.array([482.521,1226.768,1256.071])
            		authC_P2 = np.array([490.189,1087.538,746.265])
            		authC_P3 = np.array([438.284,995.631,707.942])
            		ref = np.array([authC_P0, authC_P1, authC_P2, authC_P3])
            
            if (nup == "NUP43"):
            	if (term == "N"):
            		authN_Q0 = np.array([566.331,1160.880,1273.265])
            		authN_Q1 = np.array([548.540,1269.649,1297.437])
            		authN_Q2 = np.array([541.148,1023.058,710.199])
            		authN_Q3 = np.array([481.398,931.049,667.126])
            		ref = np.array([authN_Q0, authN_Q1, authN_Q2, authN_Q3])
            
            	if (term == "C"):
            		authC_Q0 = np.array([562.660,1164.450,1271.629])
            		authC_Q1 = np.array([545.559,1273.324,1294.617])
            		authC_Q2 = np.array([536.415,1021.024,711.899])
            		authC_Q3 = np.array([476.692,928.727,668.988])
            		ref = np.array([authC_Q0, authC_Q1, authC_Q2, authC_Q3])
            
            if (nup == "NUP160"):
            	if (term == "N"):
            		authN_R0 = np.array([431.417,1128.924,1187.625])
            		authN_R1 = np.array([419.858,1238.637,1206.937])
            		authN_R2 = np.array([419.063,1089.579,783.878])
            		authN_R3 = np.array([374.947,1007.688,757.508])
            		ref = np.array([authN_R0, authN_R1, authN_R2, authN_R3])
            
            	if (term == "C"):
            		authC_R0 = np.array([530.835,1057.266,1303.285])
            		authC_R1 = np.array([503.042,1170.174,1336.088])
            		authC_R2 = np.array([554.330,1117.970,680.674])
            		authC_R3 = np.array([479.757,1041.672,632.689])
            		ref = np.array([authC_R0, authC_R1, authC_R2, authC_R3])
            
            if (nup == "NUP37"):
            	if (term == "N"):
            		authN_S0 = np.array([499.551,1124.058,1203.794])
            		authN_S1 = np.array([482.126,1240.400,1238.373])
            		authN_S2 = np.array([484.475,1085.327,770.158])
            		authN_S3 = np.array([430.940,980.348,727.628])
            		ref = np.array([authN_S0, authN_S1, authN_S2, authN_S3])
            
            	if (term == "C"):
            		authC_S0 = np.array([504.133,1142.462,1205.194])
            		authC_S1 = np.array([484.010,1259.417,1239.955])
            		authC_S2 = np.array([486.844,1066.670,772.726])
            		authC_S3 = np.array([427.069,961.708,730.504])
            		ref = np.array([authC_S0, authC_S1, authC_S2, authC_S3])
            
            if (nup == "ELYS"):
            	if (term == "N"):
            		authN_T0 = np.array([407.921,1032.623,797.957])
            		authN_T1 = np.array([348.671,956.934,773.247])
            		ref = np.array([authN_T0, authN_T1])
            
            	if (term == "C"):
            		authC_T0 = np.array([474.510,1100.386,702.278])
            		authC_T1 = np.array([412.721,1001.448,663.286])
            		ref = np.array([authC_T0, authC_T1])
            
            if (nup == "NUP98"):
            	if (term == "N"):
            		authN_U0 = np.array([646.305,1132.428,1177.901])
            		authN_U1 = np.array([514.306,1113.445,1175.581])
            		authN_U2 = np.array([631.144,1143.047,1019.882])
            		authN_U3 = np.array([612.191,1067.478,1068.629])
            		authN_U4 = np.array([605.737,1048.483,930.916])
            		authN_U5 = np.array([628.342,1121.209,872.500])
            		authN_U6 = np.array([506.732,1094.275,790.322])
            		ref = np.array([authN_U0, authN_U1, authN_U2, authN_U3, authN_U4, authN_U5, authN_U6])
            
            	if (term == "C"):
            		authC_U0 = np.array([644.062,1129.159,1204.458])
            		authC_U1 = np.array([504.275,1133.801,1157.784])
            		authC_U2 = np.array([630.039,1169.452,1031.434])
            		authC_U3 = np.array([604.717,1087.693,1087.795])
            		authC_U4 = np.array([591.400,1027.944,916.627])
            		authC_U5 = np.array([608.633,1109.125,855.262])
            		authC_U6 = np.array([488.211,1089.386,811.887])
            		ref = np.array([authC_U0, authC_U1, authC_U2, authC_U3, authC_U4, authC_U5, authC_U6])
            
            if (nup == "NUP214"):
            	if (term == "N"):
            		authN_V0 = np.array([664.617,1121.406,1258.634])
            		ref = np.array([authN_V0])
            
            	if (term == "C"):
            		authC_V0 = np.array([664.071,1135.196,1260.089])
            		ref = np.array([authC_V0])
            
            if (nup == "NUP88"):
            	if (term == "N"):
            		authN_W0 = np.array([620.590,1145.165,1267.266])
            		ref = np.array([authN_W0])
            
            	if (term == "C"):
            		authC_W0 = np.array([616.087,1200.409,1267.968])
            		ref = np.array([authC_W0])
            
                        
                        
                        
                        
            
            
        if (model == "5A9Q"):
            if (nup == "NUP43"):
                if (term == "N"):
                    authN_0 = np.array([1196.307, 578.280, 1007.303])
                    authN_R = np.array([1147.239, 748.554, 456.064])
                    authN_9 = np.array([1152.966, 676.288, 983.055])
                    authN_I = np.array([1188.192, 850.048, 420.6544])
                    ref = np.array([authN_0, authN_R, authN_9, authN_I])     
             # Nup43 C getcr #1.1/0,R,9,I:380@CA
                if (term == "C"):            
                    authC_0 = np.array([1186.279, 583.369, 1005.809])
                    authC_R = np.array([1137.054, 744.396, 458.835])
                    authC_9 = np.array([1142.872, 680.005, 979.448])
                    authC_I = np.array([1178.411, 844.515, 422.209])
                    ref = np.array([authC_0, authC_R, authC_9, authC_I])
        
            if (nup == "NUP160"):
                if (term == "N"):
            # Nup160 N getcr #1.1/1,S,a,J:41@CA
                    authN_1 = np.array([1288.803, 610.106, 903.284])
                    authN_S = np.array([1252.870, 725.826, 537.068])
                    authN_a = np.array([1260.478, 703.239, 891.066])
                    authN_J = np.array([1292.376, 823.203, 513.853])
                    ref = np.array([authN_1, authN_S, authN_a, authN_J])
                if (term == "C"):
                    # Nup160 C getcr #1.1/1,S,a,J:1201@CA
                    authC_1 = np.array([1191.789, 611.125, 983.108])
                    authC_S = np.array([1138.511, 708.726, 470.174])
                    authC_a = np.array([1146.753, 724.288, 960.323])
                    authC_J = np.array([1187.019, 805.358, 443.700])
                    ref = np.array([authC_1, authC_S, authC_a, authC_J])
    
    
            if (nup == "NUP37"):
                if(term == "N"):
                    # Nup37 N getcr #1.1/2,T,b,K:9@CA
                    authN_2 = np.array([1237.157, 598.816, 976.333])
                    authN_T = np.array([1179.777, 724.079, 484.290])
                    authN_b = np.array([1185.984, 706.003, 941.802])
                    authN_K = np.array([1226.971, 822.146, 451.789])
                    ref = np.array([authN_2, authN_T, authN_b, authN_K])
                    
                if(term == "C"):
                    # Nup37 C getcr #1.1/2,T,b,K:326@CA
                    authC_2 = np.array([1236.402, 588.357, 965.286])
                    authC_T = np.array([1181.011, 735.291, 494.523])
                    authC_b = np.array([1187.724, 694.322, 932.184])
                    authC_K = np.array([1225.765, 831.973, 463.361])
                    ref = np.array([authC_2, authC_T, authC_b, authC_K])
        
            if(nup == "SEC13"):
                if(term == "N"):
              #Sec13 N getcr #1.1/6,X,F,O:14@CA   
                    authN_6 = np.array([1140.803, 763.184, 987.984])
                    authN_X = np.array([1197.434, 664.774, 1018.729])
                    authN_F = np.array([1192.476, 762.622, 409.475])
                    authN_O = np.array([1138.455, 663.969, 443.912])
                    ref = np.array([authN_6, authN_X, authN_F, authN_O])               
                if(term == "C"):
                    authC_6 = np.array([1135.524, 781.493, 987.406])
                    authC_X = np.array([1193.033, 683.319, 1018.399])
                    authC_F = np.array([1188.354, 744.028, 410.265])
                    authC_O = np.array([1133.460, 645.602, 444.980])
                    ref = np.array([authC_6, authC_X, authC_F, authC_O])
        #Sec13 C getcr #1.1/6,X,F,O:304@CA
    
            if(nup == "SEH1"):
                if(term == "N"):
        #Seh1 N getcr #1.1/7,Y,G,P:1@CA
                    authN_7 = np.array([1135.180, 687.944, 1024.433])
                    authN_Y = np.array([1188.146, 589.668, 1051.665])
                    authN_G = np.array([1181.088, 836.930, 376.589])
                    authN_P = np.array([1131.479, 740.769, 412.990])
                    ref = np.array([authN_7, authN_Y, authN_G, authN_P])
                if(term == "C"):
        # Seh1 C getcr #1.1/7,Y,G,P:322@CA
                    authC_7 = np.array([1116.284, 703.552, 1035.049])
                    authC_Y = np.array([1173.701, 606.808, 1066.186])
                    authC_G = np.array([1167.666, 818.661, 362.467])
                    authC_P = np.array([1113.311, 726.060, 400.072])
                    ref = np.array([authC_7, authC_Y, authC_G, authC_P])
        
            if(nup == "NUP85"):
                if(term == "N"):   # Nup85 N getcr #1.1/Q,8,Z,H:9@CA
                    authN_Q = np.array([1104.998, 733.983, 412.939])
                    authN_8 = np.array([1108.652, 694.300, 1022.656])
                    authN_Z = np.array([1162.782, 599.234, 1055.191])
                    authN_H = np.array([1156.268, 826.013, 373.120])
                    ref = np.array([authN_Q, authN_8, authN_Z, authN_H])
    
                        
                if(term == "C"): # Nup85 C Atom getcrd #1.1/Q,8,Z,H:475@CA
                    authC_Q = np.array([1139.629, 717.241, 431.447])
                    authC_8 = np.array([1143.905, 709.705, 1004.162])
                    authC_Z = np.array([1195.665, 611.320, 1030.803])
                    authC_H = np.array([1189.437, 816.319, 398.179])
                    ref = np.array([authC_Q, authC_8, authC_Z, authC_H])
    
    # Nup155 N getcr #1.1/A,B:20@CA
            if(nup == "NUP155"):
                if(term == "N"):
                    authN_A = np.array([1126.810, 704.682, 588.647])
                    authN_B = np.array([1125.768, 710.075, 833.438])
                    ref = np.array([authN_A, authN_B])
    
    # Nup 155 C getcr #1.1/A,B:730@CA
                if(term == "C"):
                    authC_A = np.array([1109.824, 708.842, 587.688])
                    authC_B = np.array([1108.682, 706.273, 834.057])
                    ref = np.array([authC_A, authC_B])
               
            if nup == "NUP133": # Most N-terminal
                if (term == "N"):
                    auth3 = np.array([1145.276, 1042.219, 915.997]) #Atom #1.1/3:518@CA 
                    authU = np.array([1243.443, 946.836, 960.514]) #Atom #1.1/U:518@CA 
                    authC = np.array([1238.959, 480.349, 470.415]) #Atom #1.1/C:518@CA 
                    authL = np.array([1139.371, 391.992, 512.747]) #Atom #1.1/L:518@CA 
                    ref = np.array([auth3, authU, authC, authL]) # angles might be wrong
                
                if(term == "C"):    # Nup133 C getcr #1.1/3,U,C,L:1156@CA
                    authC_3 = np.array([1202.376, 949.301, 952.715])
                    authC_U = np.array([1270.591, 837.073, 984.803])
                    authC_C = np.array([1265.747, 584.752, 449.213])
                    authC_L = np.array([1201.188, 478.470, 477.475])
                    ref = np.array([authC_3, authC_U, authC_C, authC_L])
                    
            
            if nup == "NUP96":
                if(term == "N"): # N terminus getcr #1.1/5,W,E,N:277@CA 
                    authN_5 = np.array([1140.382, 748.182, 998.210])
                    authN_W = np.array([1196.973, 649.829, 1029.073])
                    authN_E = np.array([1191.410, 777.390, 398.923])
                    authN_N = np.array([1137.431, 678.795, 433.476])
                    ref = np.array([authN_5, authN_W, authN_E, authN_N])
    
                if(term == "C"):#Nup96 getcr getcr #1.1/5,W,E,N:751@CA 
                    auth5 = np.array([1177.224, 782.184, 975.888])    
                    authE = np.array([1229.633, 746.212, 422.979])
                    authW = np.array([1233.821, 681.917, 1004.085])
                    authN = np.array([1175.575, 645.661, 454.899])
                    ref = np.array([auth5, authE, authW, authN])
            
                
            #"select #1.1-end/M,D,V,4:150"
            if nup == "NUP107": 
                if(term == "N"):
                    authM = np.array([1134.432, 521.102, 440.365])
                    authD = np.array([1195.877, 619.745, 406.082])
                    authV = np.array([1199.033, 807.615, 1024.308])
                    auth4 = np.array([1134.991, 905.906, 994.000])
                    ref = np.array([authM, authD, authV, auth4]) 
                
                if(term == "C"): #getcr #1.1/M,D,V,4:924@CA
                    authC_4 = np.array([1162.799, 987.654, 965.130])
                    authC_V = np.array([1243.897, 883.264, 1003.383])
                    authC_D = np.array([1235.108, 543.326, 426.048])
                    authC_M = np.array([1162.346, 441.103, 460.547])
                    ref = np.array([authC_4, authC_V, authC_D, authC_M])
            
            
        # MODEL 7PEQ
        if model == "7PEQ":
            if nup == "NUP133":
                #getcr /?C:518@CA
                if(term == "N"):
                    authN_AC = np.array([541.999, 1346.763, 823.525])
                    authN_BC = np.array([638.117, 1440.660, 862.794])
                    authN_CC = np.array([523.924, 970.110, 1319.985])
                    authN_DC = np.array([591.148, 851.023, 1279.022])
                    ref = np.array([authN_AC, authN_BC, authN_CC, authN_DC])
                    
                if term == "C":
                    #getcrd /?C:1156@CA
                    authC_AC = np.array([519.158, 1237.688, 796.818])
                    authC_BC = np.array([594.420, 1338.527, 834.656])
                    authC_CC = np.array([521.047, 1080.295, 1351.747])
                    authC_DC = np.array([567.701, 958.772, 1310.625])
                    ref = np.array([authC_AC, authC_BC, authC_CC, authC_DC])
                    
            if nup == "NUP107":
                #getcrd /?D:150@CA
                if(term == "N"):
                    authN_AD = np.array([600.839, 1200.739, 781.033])
                    authN_BD = np.array([666.899, 1288.639, 811.276])
                    authN_CD = np.array([606.785, 1099.072, 1375.937])
                    authN_DD = np.array([645.965, 989.664, 1345.110])
                    ref = np.array([authN_AD, authN_BD, authN_CD, authN_DD])
                 
                #getcrd /?D:924@CAAtom 
                if(term == "C"):
                    authC_AD = np.array([551.473, 1282.074, 784.417])
                    authC_BD = np.array([633.087, 1377.111, 821.083])
                    authC_CD = np.array([542.705, 1029.732, 1363.691])
                    authC_DD = np.array([595.238, 911.805, 1324.919])
                    ref = np.array([authC_AD, authC_BD, authC_CD, authC_DD])
                    
            if nup == "NUP96":
                if(term == "N"):
                    #getcrd /?E:333@CA
                    authN_AE = np.array([581.129, 1089.071, 800.177])
                    authN_BE = np.array([631.393, 1180.644, 831.500])
                    authN_CE = np.array([620.007, 1211.549, 1351.841])
                    authN_DE = np.array([640.068, 1104.186, 1324.706])
                    ref = np.array([authN_AE, authN_BE, authN_CE, authN_DE])
                if term == "C":
                    #getcrd /?E:922@CA
                    authC_AE = np.array([556.616, 1016.471, 797.658])
                    authC_BE = np.array([592.129, 1113.714, 829.995])
                    authC_CE = np.array([596.265, 1286.142, 1356.496])
                    authC_DE = np.array([619.069, 1179.706, 1320.086])
                    ref = np.array([authC_AE, authC_BE, authC_CE, authC_DE])
                    
            if nup == "SEC13":
                if(term == "N"):
                    #getcrd /?F:14@CA
                    authN_AF = np.array([600.079, 1050.919, 795.351])
                    authN_BF = np.array([643.468, 1140.012, 825.169])
                    authN_CF = np.array([640.796, 1248.753, 1356.790])
                    authN_DF = np.array([662.173, 1140.940, 1325.852])
                    ref = np.array([authN_AF, authN_BF, authN_CF, authN_DF])
                if term == "C":
                    #getcrd /?F:304@CA
                    authC_AF = np.array([601.219, 1069.828, 795.637])
                    authC_BF = np.array([647.763, 1158.457, 825.743])
                    authC_CF = np.array([641.149, 1229.836, 1355.791])
                    authC_DF = np.array([661.812, 1122.034, 1327.058])
                    ref = np.array([authC_AF, authC_BF, authC_CF, authC_DF])
                    
            if nup == "SEH1":
                if term == "N":
                    #getcrd /?G:5@CA
                    authN_AG = np.array([627.178, 991.967, 749.919])
                    authN_BG = np.array([658.435, 1078.524, 777.520])
                    authN_CG = np.array([673.007, 1304.920, 1402.417])
                    authN_DG = np.array([695.554, 1200.844, 1365.491])
                    ref = np.array([authN_AG, authN_BG, authN_CG, authN_DG])
                
                if term == "C":
                    #getcrd /?G:320@CA
                    authC_AG = np.array([638.093, 1000.623, 749.846])
                    authC_BG = np.array([670.636, 1085.255, 777.118])
                    authC_CG = np.array([683.546, 1295.839, 1401.506])
                    authC_DG = np.array([705.751, 1191.343, 1365.910])
                    ref = np.array([authC_AG, authC_BG, authC_CG, authC_DG])     
                    
            if nup == "NUP85":
                if term == "N":
                    #getcrd /?H:20@CA
                    authN_AH = np.array([636.843, 1010.590, 740.489])
                    authN_BH = np.array([670.745, 1095.561, 768.025])
                    authN_CH = np.array([682.540, 1285.526, 1410.533])
                    authN_DH = np.array([704.151, 1182.196, 1376.040])
                    ref = np.array([authN_AH, authN_BH, authN_CH, authN_DH])
                    
                if term == "C":
                    #getcrd /?H:651@CA
                    authC_AH = np.array([570.694, 994.023, 825.221])
                    authC_BH = np.array([606.277, 1088.162, 855.288])
                    authC_CH = np.array([612.060, 1307.675, 1330.730])
                    authC_DH = np.array([636.318, 1197.414, 1292.414])
                    ref = np.array([authC_AH, authC_BH, authC_CH, authC_DH])
                    
            if nup == "NUP43":
                if term == "N":
                    #getcrd /?I:4@CA
                    authN_AI = np.array([607.400, 977.055, 801.079])
                    authN_BI = np.array([638.557, 1065.840, 829.237])
                    authN_CI = np.array([650.766, 1322.411, 1353.145])
                    authN_DI = np.array([675.070, 1213.263, 1313.972])
                    ref = np.array([authN_AI, authN_BI, authN_CI, authN_DI])
                    
                if term == "C":
                    #getcrd /?I:380@CA
                    authC_AI = np.array([615.244, 984.988, 802.829])
                    authC_BI = np.array([647.684, 1072.310, 830.776])
                    authC_CI = np.array([658.171, 1314.277, 1350.619])
                    authC_DI = np.array([682.206, 1204.628, 1312.594])
                    ref = np.array([authC_AI, authC_BI, authC_CI, authC_DI])
                    
            if nup == "NUP160":
                if term == "N":
                    #getcrd /?J:78@CA
                    authN_AJ = np.array([521.159, 984.775, 890.577])
                    authN_BJ = np.array([558.559, 1085.539, 922.612])
                    authN_CJ = np.array([558.979, 1321.315, 1268.981])
                    authN_DJ = np.array([545.739, 1234.093, 1248.102])
                    ref = np.array([authN_AJ, authN_BJ, authN_CJ, authN_DJ])
                    
                if term == "C":
                    #getcrd /?J:1195@CA
                    authC_AJ = np.array([590.711, 1029.030, 807.313])
                    authC_BJ = np.array([630.923, 1119.928, 837.190])
                    authC_CJ = np.array([632.207, 1271.078, 1346.311])
                    authC_DJ = np.array([654.672, 1162.184, 1312.629])
                    ref = np.array([authC_AJ, authC_BJ, authC_CJ, authC_DJ])
                    
            if nup == "NUP37":
                if term == "N":
                    #getcrd /?K:18@CA
                    authN_AK = np.array([533.337, 992.853, 833.045])
                    authN_BK = np.array([569.551, 1092.930, 864.748])
                    authN_CK = np.array([574.306, 1310.683, 1325.339])
                    authN_DK = np.array([590.786, 1213.351, 1280.732])
                    ref = np.array([authN_AK, authN_BK, authN_CK, authN_DK])
                    
                if term == "C":
                    #getcrd /?K:324@CA
                    authC_AK = np.array([553.058, 985.968, 839.464])
                    authC_BK = np.array([588.090, 1082.696, 870.158])
                    authC_CK = np.array([593.860, 1317.043, 1317.931])
                    authC_DK = np.array([600.402, 1212.340, 1261.130])
                    ref = np.array([authC_AK, authC_BK, authC_CK, authC_DK])
                    

        if model == "7PER":                    
            if nup == "NUP205":
                if term == "N":
                    #getcrd /D,J,V,P:9@CA
                    authN_D = np.array([766.579, 1044.574, 1052.258])
                    authN_J = np.array([708.820, 1125.978, 1045.733])
                    authN_V = np.array([686.322, 1092.394, 1186.438])
                    authN_P = np.array([753.571, 1153.951, 1195.119])
                    ref = np.array([authN_D, authN_J, authN_V, authN_P])
                                    
                if term == "C":       
                    #getcrd /D,J,V,P:1692@CA
                    authC_D = np.array([683.834, 1148.509, 1097.252])
                    authC_J = np.array([742.309, 1231.826, 1045.529])
                    authC_V = np.array([717.704, 986.583, 1198.460])
                    authC_P = np.array([670.166, 1058.951, 1134.355])
                    ref = np.array([authC_D, authC_J, authC_V, authC_P])
                    
            if nup == "NUP155":
                if term == "N":
                    #getcrd /E,K,Q,W:20@CA
                    authN_E = np.array([668.125, 986.789, 1118.243])
                    authN_K = np.array([649.617, 1074.233, 1045.931])
                    authN_Q = np.array([676.362, 1246.211, 1104.215])
                    authN_W = np.array([638.337, 1131.825, 1177.548])
                    ref = np.array([authN_E, authN_K, authN_Q, authN_W])
                if term == "C":                    
                    #getcrd /E,K,Q,W:1375@CA
                    authC_E = np.array([701.274, 1086.404, 1084.926])
                    authC_K = np.array([718.793, 1165.954, 1023.777])
                    authC_Q = np.array([708.040, 1150.287, 1154.007])
                    authC_W = np.array([698.384, 1045.048, 1205.549])
                    ref = np.array([authC_E, authC_K, authC_Q, authC_W])
                    
            if nup == "NUP93":
                if term == "N":
                    #getcrd /C,I,O,U:173@CA
                    authN_C = np.array([679.385, 1115.497, 1107.962])
                    authN_I = np.array([693.222, 1207.092, 1063.216])
                    authN_O = np.array([683.229, 1104.768, 1124.793])
                    authN_U = np.array([672.276, 1012.543, 1162.728])
                    ref = np.array([authN_C, authN_I, authN_O, authN_U])
                if term == "C":     
                    #getcrd /C,I,O,U:815@CA
                    authC_C = np.array([664.864, 1048.350, 1112.871])
                    authC_I = np.array([663.228, 1147.690, 1045.459])
                    authC_O = np.array([669.738, 1171.469, 1114.178])
                    authC_U = np.array([633.985, 1069.780, 1161.499])
                    ref = np.array([authC_C, authC_I, authC_O, authC_U])
                    
            if nup == "NUP54":
                if term == "N":
                    #getcrd /F,L,R,X:128@CA
                    authN_F = np.array([793.780, 1098.840, 1103.565])
                    authN_L = np.array([796.577, 1199.915, 1102.008])
                    authN_R = np.array([785.953, 1096.559, 1151.106])
                    authN_X = np.array([782.070, 996.679, 1138.960])
                    ref = np.array([authN_F, authN_L, authN_R, authN_X])
                    
                if term == "C":   
                    #getcrd /F,L,R,X:493@CA
                    authC_F = np.array([733.940, 1032.037, 1143.314])
                    authC_L = np.array([744.925, 1117.484, 1114.681])
                    authC_R = np.array([741.041, 1166.010, 1098.352])
                    authC_X = np.array([746.295, 1087.568, 1129.866])
                    ref = np.array([authC_F, authC_L, authC_R, authC_X])
                    
            if nup == "NUP58":
                if term == "N":
                    #getcrd /G,M,S,Y:248@CA
                    authN_G = np.array([796.774, 1096.951, 1088.269])
                    authN_M = np.array([805.350, 1200.433, 1088.999])
                    authN_S = np.array([787.615, 1099.730, 1166.393])
                    authN_Y = np.array([789.007, 993.510, 1152.683])
                    ref = np.array([authN_G, authN_M, authN_S, authN_Y])
                if term == "C":    
                    #getcrd /G,M,S,Y:418@CA
                    authC_G = np.array([741.218, 1036.293, 1147.267])
                    authC_M = np.array([748.939, 1122.796, 1121.192])
                    authC_S = np.array([747.976, 1160.386, 1095.709])
                    authC_Y = np.array([750.157, 1082.239, 1123.280])
                    ref = np.array([authC_G, authC_M, authC_S, authC_Y])
                    
            if nup == "NUP62":
                if term == "N":
                    #getcrd /H,N,T,Z:334@CA
                    authN_H = np.array([782.857, 1100.463, 1096.103])
                    authN_N = np.array([789.127, 1199.869, 1090.957])
                    authN_T = np.array([774.205, 1097.290, 1157.357])
                    authN_Z = np.array([773.444, 997.033, 1149.114])
                    ref = np.array([authN_H, authN_N, authN_T, authN_Z])
                if term == "C":       
                    #getcrd /H,N,T,Z:502@CA
                    authC_H = np.array([727.667, 1037.177, 1144.194])
                    authC_N = np.array([737.704, 1120.965, 1113.178])
                    authC_T = np.array([734.202, 1161.708, 1097.243])
                    authC_Z = np.array([738.479, 1085.264, 1130.228])
                    ref = np.array([authC_H, authC_N, authC_T, authC_Z])
                    
        if model == "5IJN":
            if nup == "NUP205":
                if term == "N":
                    #getcrd #1.1/D,J,P,V:9@CA
                    authN_D = np.array([960.284, 652.509, 782.339])
                    authN_J = np.array([961.280, 678.281, 811.489])
                    authN_P = np.array([960.962, 767.045, 644.009])
                    authN_V = np.array([964.721, 743.088, 614.105])
                    ref = np.array([authN_D, authN_J, authN_P, authN_V])    
                    
                if term == "C":
                    #getcrd #1.1/D,J,P,V:1692@CA
                    authC_D = np.array([1046.501, 746.245, 721.493])
                    authC_J = np.array([992.239, 816.461, 788.789])
                    authC_P = np.array([1038.535, 667.639, 704.642])
                    authC_V = np.array([990.301, 598.984, 634.525])
                    ref = np.array([authC_D, authC_J, authC_P, authC_V])
                    
            if nup == "NUP54":
                if term == "N":
                    #getcrd #1.1/F,L,R,X:128@CA
                    authN_F = np.array([934.656, 716.059, 736.701])
                    authN_L = np.array([938.484, 806.287, 742.363])
                    authN_R = np.array([932.467, 700.520, 695.076])
                    authN_X = np.array([935.376, 612.148, 686.866])
                    ref = np.array([authN_F, authN_L, authN_R, authN_X])   
                    
                if term == "C":
                    #getcrd #1.1/F,L,R,X:493@CA
                    authC_F = np.array([984.530, 653.126, 680.352])
                    authC_L = np.array([970.364, 714.857, 726.627])
                    authC_R = np.array([987.751, 762.904, 746.799])
                    authC_X = np.array([967.570, 703.607, 701.773])
                    ref = np.array([authC_F, authC_L, authC_R, authC_X])

                
            if nup == "NUP58":
                if term == "N":
                    #getcrd #1.1/G,M,S,Y:248@CA
                    authN_G = np.array([930.053, 710.449, 750.623])
                    authN_M = np.array([930.742, 808.497, 755.842])
                    authN_S = np.array([927.096, 706.032, 681.392])
                    authN_Y = np.array([928.412, 609.629, 673.023])
                    ref = np.array([authN_G, authN_M, authN_S, authN_Y])
                if term == "C":     
                    #getcrd #1.1/G,M,S,Y:418@CA
                    authC_G = np.array([978.520, 659.710, 677.663])
                    authC_M = np.array([967.247, 721.001, 720.364])
                    authC_S = np.array([981.755, 756.542, 750.005])
                    authC_Y = np.array([964.003, 697.693, 708.018])
                    ref = np.array([authC_G, authC_M, authC_S, authC_Y])
            
            if nup == "NUP155":
                if term == "N":
                    #getcrd #1.1/A,B,E,K,Q,W:20@CA
                    authN_A = np.array([1127.765, 714.280, 844.037])
                    authN_B = np.array([1129.475, 701.251, 582.893])
                    authN_E = np.array([1035.081, 602.042, 702.602])
                    authN_K = np.array([1091.232, 679.306, 777.270])
                    authN_Q = np.array([1037.108, 812.567, 723.198])
                    authN_W = np.array([1093.999, 718.479, 650.195])
                    ref = np.array([authN_A, authN_B, authN_E, authN_K, authN_Q, authN_W])
                    
                if term == "C": 
                    #getcrd #1.1/A,B:863@CA, getcrd #1.1/E,K,Q,W:1375@CA 
                    authC_A = np.array([1123.935, 733.451, 853.606])
                    authC_B = np.array([1126.827, 681.928, 573.229])
                    authC_E = np.array([1050.620, 694.556, 764.925])
                    authC_K = np.array([1033.018, 765.870, 830.250])
                    authC_Q = np.array([1051.867, 722.745, 661.181])
                    authC_W = np.array([1031.677, 647.801, 595.052])
                    ref = np.array([authC_A, authC_B, authC_E, authC_K, authC_Q, authC_W])
            
            if nup == "NUP93":
                if term == "N":
                    #getcrd #1.1/C,I,O,U:1@CA
                    authN_C = np.array([976.397, 666.536, 729.429])
                    authN_I = np.array([961.761, 749.355, 763.917])
                    authN_O = np.array([975.962, 748.926, 698.630])
                    authN_U = np.array([960.689, 668.109, 665.067])
                    ref = np.array([authN_C, authN_I, authN_O, authN_U])
                    
                if term == "C": 
                    #getcrd #1.1/C,I,O,U:815@CA
                    authC_C = np.array([1060.023, 674.793, 712.464])
                    authC_I = np.array([1055.745, 758.808, 794.813])
                    authC_O = np.array([1057.252, 740.183, 714.214])
                    authC_U = np.array([1054.300, 653.797, 630.882])
                    ref = np.array([authC_C, authC_I, authC_O, authC_U])

            if nup == "NUP62":
                if term == "N":
                    #getcrd #1.1/H,N,T,Z:334@CA
                    authN_H = np.array([945.004, 713.228, 744.610])
                    authN_N = np.array([946.332, 804.477, 752.982])
                    authN_T = np.array([942.333, 702.920, 686.442])
                    authN_Z = np.array([943.862, 613.538, 676.683])
                    ref = np.array([authN_H, authN_N, authN_T, authN_Z])
                    
                if term == "C": 
                    #getcrd #1.1/H,N,T,Z:502@CA
                    authC_H = np.array([991.687, 656.938, 681.239])
                    authC_N = np.array([978.236, 716.679, 727.738])
                    authC_T = np.array([994.722, 758.872, 745.506])
                    authC_Z = np.array([975.468, 701.654, 701.185])
                    ref = np.array([authC_H, authC_N, authC_T, authC_Z])

                                     
        if model == "5IJO":                    
            if nup == "NUP54":
                if term == "N":
                    #getcrd /F,L,R,X:128@CA
                    authN_F = np.array([934.656, 716.059, 736.701])
                    authN_L = np.array([938.484, 806.287, 742.363])
                    authN_R = np.array([932.467, 700.520, 695.076])
                    authN_X = np.array([935.376, 612.148, 686.866])
                    ref = np.array([authN_F, authN_L, authN_R, authN_X])

                if term == "C": 
                    #getcrd /F,L,R,X:493@CA
                    authC_F = np.array([984.530, 653.126, 680.352])
                    authC_L = np.array([970.364, 714.857, 726.627])
                    authC_R = np.array([987.751, 762.904, 746.799])
                    authC_X = np.array([967.570, 703.607, 701.773])
                    ref = np.array([authC_F, authC_L, authC_R, authC_X])

            if nup == "NUP188":
                if term == "N":
                    #getcrd /J,V:1@CA
                    authN_J = np.array([963.691, 696.108, 808.948])
                    authN_V = np.array([964.642, 722.851, 612.857])
                    ref = np.array([authN_J, authN_V])
                    
                if term == "C": 
                    #getcrd /J,V:1564@CA
                    authC_J = np.array([989.042, 805.579, 797.104])
                    authC_V = np.array([987.194, 611.748, 624.551])
                    ref = np.array([authC_J, authC_V])
                    
            if nup == "NUP205":
                if term == "N":
                    #getcrd /D,P:9@CA
                    authN_D = np.array([960.284, 652.509, 782.339])
                    authN_P = np.array([960.962, 767.045, 644.009])
                    ref = np.array([authN_D, authN_P])
                    
                if term == "C": 
                    #getcrd /D,P:1692@CA
                    authC_D = np.array([1046.501, 746.245, 721.493])
                    authC_P = np.array([1038.535, 667.639, 704.642])
                    ref = np.array([authC_D, authC_P])

            if nup == "NUP155":
                if term == "N":
                    #getcrd /A,B,E,K,Q,W:20@CA
                    authN_A = np.array([1127.765, 714.280, 844.037]) # as in 5IJN 1.1
                    authN_B = np.array([1129.475, 701.251, 582.893])
                    authN_E = np.array([1035.081, 602.042, 702.602])
                    authN_K = np.array([1091.232, 679.306, 777.270])
                    authN_Q = np.array([1037.108, 812.567, 723.198])
                    authN_W = np.array([1093.999, 718.479, 650.195])
                    ref = np.array([authN_A, authN_B, authN_E, authN_K, authN_Q, authN_W])
                    
                if term == "C": 
                    #getcrd /A,B:863@CA, getcrd /E,K,Q,W:1375@CA
                    authC_A = np.array([1123.935, 733.451, 853.606])
                    authC_B = np.array([1126.827, 681.928, 573.229])
                    authC_E = np.array([1050.620, 694.556, 764.925])
                    authC_K = np.array([1033.018, 765.870, 830.250])
                    authC_Q = np.array([1051.867, 722.745, 661.181])
                    authC_W = np.array([1031.677, 647.801, 595.052])
                    ref = np.array([authC_A, authC_B, authC_E, authC_K, authC_Q, authC_W])
                    
            if nup == "NUP93":
                if term == "N":
                    #getcrd /C,I,O,U:1@CA
                    authN_C = np.array([976.397, 666.536, 729.429])
                    authN_I = np.array([961.761, 749.355, 763.917])
                    authN_O = np.array([975.962, 748.926, 698.630])
                    authN_U = np.array([960.689, 668.109, 665.067])
                    ref = np.array([authN_C, authN_I, authN_O, authN_U])
                    
                if term == "C": 
                    #getcrd /C,I,O,U:815@CA
                    authC_C = np.array([1060.023, 674.793, 712.464])
                    authC_I = np.array([1055.745, 758.808, 794.813])
                    authC_O = np.array([1057.252, 740.183, 714.214])
                    authC_U = np.array([1054.300, 653.797, 630.882])
                    ref = np.array([authC_C, authC_I, authC_O, authC_U])
                    
            if nup == "NUP58":
                if term == "N":
                    #getcrd /G,M,S,Y:248@CA
                    authN_G = np.array([930.053, 710.449, 750.623])
                    authN_M = np.array([930.742, 808.497, 755.842])
                    authN_S = np.array([927.096, 706.032, 681.392])
                    authN_Y = np.array([928.412, 609.629, 673.023])
                    ref = np.array([authN_G, authN_M, authN_S, authN_Y])
                    
                if term == "C":     
                    #getcrd /G,M,S,Y:418@CA
                    authC_G = np.array([978.520, 659.710, 677.663])
                    authC_M = np.array([967.247, 721.001, 720.364])
                    authC_S = np.array([981.755, 756.542, 750.005])
                    authC_Y = np.array([964.003, 697.693, 708.018])
                    ref = np.array([authC_G, authC_M, authC_S, authC_Y])
                
            if nup == "NUP62":
                if term == "N":
                    #getcrd /H,N,T,Z:334@CA
                    authN_H = np.array([945.004, 713.228, 744.610])
                    authN_N = np.array([946.332, 804.477, 752.982])
                    authN_T = np.array([942.333, 702.920, 686.442])
                    authN_Z = np.array([943.862, 613.538, 676.683])
                    ref = np.array([authN_H, authN_N, authN_T, authN_Z])
                    
                if term == "C": 
                    #getcrd /H,N,T,Z:502@CA
                    authC_H = np.array([991.687, 656.938, 681.239])
                    authC_N = np.array([978.236, 716.679, 727.738])
                    authC_T = np.array([994.722, 758.872, 745.506])
                    authC_Z = np.array([975.468, 701.654, 701.185])
                    ref = np.array([authC_H, authC_N, authC_T, authC_Z])
                    
        return ref 
    
    def centre(self, model):
        if model == "5A9Q":  
            c = np.array([711.36, 711.36]) # central axis 

        if model == "SIMPLE":  
            c = np.array([711.36, 711.36]) # central axis 
            
        if model == "OTHER" or model == "SMILE": 
            c = np.array([0, 0])
            
        if model == "7PEQ" or model == "7PER": 
            c = np.array([1094.4, 1094.4]) # central axis 7PEQ determined from EM density map
            
        if model == "7R5K":
            c = np.array([970.56, 970.56]) # from EM density map, verified with Nup107 opposite corners
            
        if model == "7R5J":
            c = np.array([993.6, 993.6])    # from EM density map, verified with Nup107 opposite corners
            
        if model == "5IJN" or model == "5IJO" :
            c = np.array([716.9, 713.55])
            
        return c

    def assignRing(self, ref, nup, model):
        if model == "5A9Q" or model == "7PEQ": 
            if nup != "NUP155":
                return ["CR" if ref[i,2] > np.mean(ref[:,2]) else "NR" for i in range(len(ref))]
            else:
                return ["BRCR" if ref[i,2] > np.mean(ref[:,2]) else "BRNR" for i in range(len(ref))]
            
        if model == "SIMPLE" or model == "OTHER" or model == "SMILE":
            return ["IR" for i in range(len(ref))]
            
        if model == "7PER" or model == "5IJN" or model == "5IJO":
            if nup != "NUP155":
                return ["IR" for i in range(len(ref))]
            else:
                return ["BRCR" if ref[i,2] == np.max(ref[:,2]) else "BRNR" if ref[i,2] == np.min(ref[:,2]) else "IR" for i in range(len(ref))]
                
        if model == "7R5J" or model == "7R5K":
            ORonly = ["NUP133", "NUP107", "NUP96", "SEC13", "SEH1", "NUP85", "NUP160", "NUP37", "NUP43"] 
            IRonly = ["NUP210", "ALADIN", "NUP188", "NDL1", "NUP35", "NUP54", "NUP58"]
            CRonly = ["RANBP2", "NUP214", "NUP88"]
            
            if nup in ORonly:
                return ["CR" if ref[i,2] > np.mean(ref[:,2]) else "NR" for i in range(len(ref))]

            elif nup in IRonly: 
                return ["IR" for i in range(len(ref))]
            
            elif nup in CRonly:
                return ["CR" for i in range(len(ref))]
            
            elif nup == "ELYS":
                return ["NR" for i in range(len(ref))]
            
            elif nup == "NUP93" or nup == "NUP205" or nup == "NUP98":
                sortz = np.sort(ref[:,2])    
                return ["NR" if ref[i, 2] <= sortz[0] else "CR" if ref[i,2] >= sortz[-2] else "IR" for i in range(len(ref))]
            
            elif nup == "NUP62":
                sortz = np.sort(ref[:,2])    
                return ["CR" if ref[i,2] >= sortz[-1] else "IR" for i in range(len(ref))]
            
            elif nup == "NUP155":
                sortz = np.sort(ref[:,2])  
                return ["BRNR" if ref[i, 2] == sortz[0] else "BRCR" if ref[i, 2] == sortz[-1] else "IR" for i in range(len(ref))]

            
            
    
    def cornang(self, p_ref, p1): #TODO: Does this also work for rotational units in a different quadrant? 
        "Angle between node p0, p1 (PDB 5A9Q, coordinates from chimera) and the centre"
        c = self.centre(self.model)
            
        p_refnew = p_ref[:2] - c[:2]   # change coordinate system to 0 centre
        p1new = p1[:2] - c[:2]
        
        rotateby = np.arctan2(p_refnew[1], p_refnew[0]) # rotate such that ref overlays with positive x axis 
        arctan2_p1 = np.arctan2(p1new[1], p1new[0])
        ang = arctan2_p1 - rotateby

        if ang < -math.pi:
            ang += 2*math.pi
        if ang > math.pi:
            ang -= 2*math.pi

#  #      if (abs(ang) > math.pi): ang = 2*math.pi - ang

        return ang#np.arctan2(p1new[1], p1new[0]) - rotateby

    
    def allang(self, ref):

        ringAngles = []
        for i in range(len(ref)):
            ringAngles.append(self.cornang(ref[0], ref[i]))
            
        ringAngles = ringAngles - np.mean(ringAngles) 
        
        return ringAngles # note: If the rotation is counterclockwise, the angle has a positive measure
    
    def distz(self, p):
        
        z = np.zeros(len(p))
    
        for i in range(len(z)):
            z[i] = p[i][2]
    
        z -= min(z) # offset so that smallest z is 0 TODO: this needs to be global for multicolour 
        z /= 10 # Angstrom to nm 
        z = [round(i, ndigits = 2) for i in z] # round to nearest .1 nm
        return np.array(z)
    
    def radii(self, p):
        c = self.centre(self.model)

        r = np.zeros(len(p))
        for i in range(len(r)):
            x = p[i][0] - c[0]
            y = p[i][1] - c[1]
            r[i] = np.sqrt(x**2 + y**2)
    
            
        r /= 10
        r = [round(i, ndigits=1) for i in r]
        return r
    
    def rotang(self, z, ringAng):
        "Find azimuthal angle between nearest 'corners' of two rings"
        midplane = np.mean(z) 

        octOffset = 2*np.pi/8 # TODO: different sym?

        crAng = np.mean(ringAng[z > midplane]) 
        crAngCW = crAng - octOffset # next rot unit CW
        crAngACW = crAng + octOffset # next rot unit ACW
        
        #Angle of rotational unit in NR
        nrAng = np.mean(ringAng[z < midplane])

        ang = [(nrAng - i) for i in [crAngCW, crAng, crAngACW]]

        minval, minindex = min([(abs(val), idx) for (idx, val) in enumerate(ang)])
        
        offset = [-octOffset, 0, octOffset][minindex]
        return ang[minindex], offset
        
        
