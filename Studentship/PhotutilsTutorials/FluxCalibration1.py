

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.io import fits
from stackMultiple import imageStack
from FinalAperturePhotometry import* 
from psfForStackedImages import*
from astropy.stats import sigma_clip, sigma_clipped_stats
import os

def OutlierRegress(xlist, ylist):

    """
    double layered regression with two different clipping parameters
    """

    regress = linregress(xlist, ylist)
    slope = regress.slope
    intercept = regress.intercept
    def y(x):
        return (slope*x + intercept)

    distList = []
    i = 0
    for i in range(0,len(xlist)):
        xpos = xlist[i]
        ypos = ylist[i]

        yLine = y(xpos)
        yDist = np.abs(yLine - ypos)
        distList.append(yDist)
        i += 1

    clipDist = sigma_clip(distList, sigma = 2.5, axis = 0, masked = False)
    
    indVal = np.argwhere(np.isnan(clipDist))
    newXlist = np.delete(xlist, indVal)
    newYlist = np.delete(ylist, indVal)
    
    regress2 = linregress(xlist, ylist)
    slope2 = regress2.slope
    intercept2 = regress2.intercept

    distList2 = []
    ii = 0
    for ii in range(0,len(newXlist)):
        xpos = newXlist[ii]
        ypos = newYlist[ii]

        yLine = y(xpos)
        yDist = np.abs(yLine - ypos)
        distList2.append(yDist)
        ii += 1

    clipDist2 = sigma_clip(distList2, sigma = 5.0, axis = 0, masked = False)
    indVal2 = np.argwhere(np.isnan(clipDist2))
    newXlist2 = np.delete(newXlist, indVal2)
    newYlist2 = np.delete(newYlist, indVal2)


    return [newXlist2, newYlist2]




def FluxCalibration(objectName, method, date, stackNum): 

    """
    objectName - single object name as a string, e.g. 'WD0830'
    method - 'mean' or 'median'
    date - string of dates or single date, in yyyymmdd format. e.g. '20221214,20230112,20230114...'
    stackNum - 8 times amount of dates
    """

    file1 = f'Studentship/PhotutilsTutorials/{objectName}Stars.txt'
    if method == 'median':
        file2 = f'Studentship/PhotutilsTutorials/{date}_{objectName}_MedianStackedApertureResults-B_{stackNum}.txt'
        file3 = f'Studentship/PhotutilsTutorials/{date}_{objectName}_MedianStackedApertureResults-V_{stackNum}.txt'
    else:
        file2 = f'Studentship/PhotutilsTutorials/{date}_{objectName}_MeanStackedApertureResults-B_{stackNum}.txt'
        file3 = f'Studentship/PhotutilsTutorials/{date}_{objectName}_MeanStackedApertureResults-V_{stackNum}.txt'

    with open(file2) as f4:
        lines4 = f4.readlines()

    i = 0
    ignoreList = []
    for line4 in lines4:
        var = line4.strip().strip('\n').split(' ')
        if float(var[3]) < 0:
            ignoreList.append(i)
        i += 1

    with open(file1, 'r') as f:
        lines = f.readlines()

    BFluxList = []
    BFluxErrorList = []
    VFluxList = []
    VFluxErrorList = []

    knownMags_B = []
    knownMags_V = []
    j = 0
    for line in lines:
        var = line.strip('\n').split(' ')
        if j not in ignoreList:
            knownMags_B.append(float(var[2]))
            knownMags_V.append(float(var[3]))
        j += 1

    with open(file2, 'r') as f2:
        lines2 = f2.readlines()

    for line2 in lines2:
        var2 = line2.strip('\n').split(' ')
        if float(var2[3]) > 0:
            BFluxList.append(float(var2[3]))
            BFluxErrorList.append(float(var2[4]))

    with open(file3, 'r') as f3:
        lines3 = f3.readlines()

    for line3 in lines3:
        var3 = line3.strip('\n').split(' ')
        if float(var3[3]) > 0:
            VFluxList.append(float(var3[3]))
            VFluxErrorList.append(float(var3[4]))

    standardMag_B = []
    standardMag_V = []
    

    knownMagB = knownMags_B[0]
    knownMagV = knownMags_V[0]
    knownFluxB = BFluxList[0]
    knownFluxB_error = BFluxErrorList[0]
    knownFluxV = VFluxList[0]
    knownFluxV_error = VFluxErrorList[0]
    standardMag_B_error = []
    standardMag_V_error = []

    i = 0
    for i in range(0, len(BFluxList)):
        targetFluxB = float(BFluxList[i])
        targetFluxB_error = float(BFluxErrorList[i])
        targetFluxV = float(VFluxList[i])
        targetFluxV_error = float(VFluxErrorList[i])
        standardMagTarget_B = knownMagB -2.5*np.log10(targetFluxB / knownFluxB)
        standardMagTarget_B_error = np.sqrt((-2.5 * targetFluxB_error / (np.log(10) * targetFluxB ) )**2 + ( 2.5 * knownFluxB_error / (np.log(10) * knownFluxB) )**2)
        standardMagTarget_V = knownMagV -2.5*np.log10(targetFluxV / knownFluxV)
        standardMagTarget_V_error = np.sqrt(((-2.5 * targetFluxV_error) / (np.log(10) * targetFluxV ) )**2 + ( (2.5 * knownFluxV_error) / (np.log(10) * knownFluxV) )**2)
        standardMag_B.append(standardMagTarget_B)
        standardMag_V.append(standardMagTarget_V)
        standardMag_B_error.append(standardMagTarget_B_error)
        standardMag_V_error.append(standardMagTarget_V_error)
    
    instrumentalMagB = []
    instrumentalMagV = []
    for jj in range(0,len(BFluxList)):
        instrumentalMagB.append(-2.5 * np.log10(float(BFluxList[jj] / 90)))
        instrumentalMagV.append(-2.5 * np.log10(float(VFluxList[jj] / 90)))
    instrumentalMagB = np.array(instrumentalMagB)
    instrumentalMagV = np.array(instrumentalMagV)
    zeroPointListB = np.abs(np.subtract(instrumentalMagB,standardMag_B))
    zeroPointListV = np.abs(np.subtract(instrumentalMagV,standardMag_V))

    compB = knownMags_B[0]
    compV = knownMags_V[0]
    compb = standardMag_B[0]
    compb_error = standardMag_B_error[0]
    compv = standardMag_V[0]
    compv_error = standardMag_V_error[0]

    standardMag_B = np.array(standardMag_B[1:]) 
    standardMag_V = np.array(standardMag_V[1:])
    knownMagsB = np.array(knownMags_B[1:])
    knownMagsV = np.array(knownMags_V[1:])

    BminusVMeasured = np.subtract(standardMag_B, standardMag_V)
    BminusVKnown = np.subtract(knownMagsB, knownMagsV)

    Bminusb = np.subtract(knownMagsB, standardMag_B)
    Vminusv = np.subtract(knownMagsV, standardMag_V)

    newXY = OutlierRegress(BminusVKnown, BminusVMeasured)
    linesParams = linregress(newXY[0], newXY[1])
    grad = float(linesParams.slope)
    graderror = linesParams.stderr
    intercept = linesParams.intercept

    x = np.array([-0.3, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4])
    y = grad*x + intercept


    Tbv = 1/grad
    Tbv_error = graderror / (grad**2)
    plt.plot(x,y)
    plt.scatter(BminusVKnown, BminusVMeasured)
    plt.xlabel("Known B-V")
    plt.ylabel("Calculated B-V")
    #plt.show()


    newXY2 = OutlierRegress(BminusVKnown, Bminusb)
    linesParams2 = linregress(newXY2[0], newXY2[1])
    Tb_bv = linesParams2.slope
    Tb_bv_error = linesParams2.stderr
    intercept2 = linesParams2.intercept

    x = np.array([-0.3, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4])
    y = Tb_bv*x + intercept2

    plt.plot(x,y)
    plt.scatter(BminusVKnown, Bminusb)
    plt.xlabel("Known B-V")
    plt.ylabel("Measured B-b")
    #plt.show()

    

    newXY3 = OutlierRegress(BminusVKnown, Vminusv)
    linesParams3 = linregress(newXY3[0], newXY3[1])
    Tv_bv = linesParams3.slope
    Tv_bv_error = linesParams3.stderr
    intercept3 = linesParams3.intercept

    x = np.array([-0.3, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4])
    y = Tv_bv*x + intercept3

    plt.plot(x,y)
    plt.scatter(BminusVKnown, Vminusv)
    plt.xlabel("Known B-V")
    plt.ylabel("Measured V-v")
    #plt.show()

    print("The value of Tbv is: ", Tbv, " +- ", Tbv_error)
    print("The value of Tb_bv is: ", Tb_bv, " +- ", Tb_bv_error)
    print("The value of Tv_bv is: ", Tv_bv, " +- ", Tv_bv_error)


    correctedMagB1 = []
    correctedMagV1 = []
    correctedMagB1_error = []
    correctedMagV1_error = []
    correctedMagB2 = []
    correctedMagV2 = []

    for ii in range(0, len(standardMag_B)):
        B_target1 = standardMag_B[ii] + compB - compb + Tb_bv * Tbv * ((standardMag_B[ii] - standardMag_V[ii]) - (compb - compv))
        B_target1_error = np.sqrt( ((Tb_bv * Tbv + 1) * standardMag_B_error[ii])**2 + ((-1*Tb_bv*Tbv - 1) * compb_error)**2 + (Tb_bv*(compv-compb + standardMag_B[ii] - standardMag_V[ii]) * Tbv_error)**2 + (Tbv*(compv-compb + standardMag_B[ii] - standardMag_V[ii]) * Tb_bv_error )**2 + (-1*Tb_bv*Tbv*standardMag_V_error[ii])**2 + (Tb_bv * Tbv * compv_error)**2 )
        V_target1 = B_target1 + Tbv * ((compb - compv)-(standardMag_B[ii] - standardMag_V[ii])) - (compB - compV)
        V_target1_error = np.sqrt( (B_target1_error)**2 + ((compb-compv-standardMag_B[ii] + standardMag_V[ii])*Tbv_error)**2 + (Tbv * compb_error)**2 + (-1*Tbv*compv_error)**2 + (-1*Tbv*standardMag_B_error[ii])**2 + (Tbv*standardMag_V_error[ii])**2)
        correctedMagB1.append(B_target1)
        correctedMagB1_error.append(B_target1_error)
        correctedMagV1.append(V_target1)
        correctedMagV1_error.append(V_target1_error)

        V_target2 = standardMag_V[ii] + (compV - compv) + Tbv * Tv_bv * ((standardMag_B[ii] - standardMag_V[ii]) - (compb - compv))
        B_target2 = V_target2 + Tbv * ((standardMag_B[ii] - standardMag_V[ii]) - (compb - compv)) + compB - compV
        correctedMagV2.append(V_target2)
        correctedMagB2.append(B_target2)

    diffB1 = np.subtract(knownMagsB, correctedMagB1)
    diffV1 = np.subtract(knownMagsV, correctedMagV1)
    Diff = [diffB1,diffV1]
    print('Known magnitudes subtract the calculated B magnitudes: ')
    print(diffB1)
    print('\nKnown magnitudes subtract the calculated V magnitudes: ')
    print(diffV1)

    diffB2 = np.subtract(knownMagsB, correctedMagB2)
    diffV2 = np.subtract(knownMagsV, correctedMagV2)

    f2 = open(f'Studentship/PhotutilsTutorials/{date}_{objectName}_{stackNum}_{method}-Coefficients.txt', 'w')
    f2.write(f'Tbv, {Tbv}, {Tbv_error}\n')
    f2.write(f'Tb_bv, {Tb_bv}, {Tb_bv_error}\n')
    f2.write(f'Tv_bv, {Tv_bv}, {Tv_bv_error}\n')
    f2.write(f'ZeroPointB, {zeroPointListB[0]}\n')
    f2.write(f'ZeroPointV, {zeroPointListV[0]}\n')
    f2.close()

    f3 = open(f'Studentship/PhotutilsTutorials/{date}_{objectName}_{stackNum}_{method}-Magnitudes.txt', 'w')
    for jj in range(0, len(correctedMagB1)):
        f3.write(f'{knownMagsB[jj]} {correctedMagB1[jj]} {correctedMagB1_error[jj]}\n')

    f3.close()

    return Diff

def RepeatedFluxCalibration(objects:str, dates:str, method, mode):
    """
    
    for objects, write in the format: object1,object2,object3,object4...
    dates would be corresponding dates: 20221214,20230114...

    METHOD is either 'mean' or 'median' will mostly be using mean

    MODE is either 'all' or 'single'. This controls whether every object/data pair is a unique 
    object, or whether they are treated as one large dataset.

    """
    if mode == 'single' or mode == 'all':
        dateList = dates.split(',')
        ObjectList = objects.split(',')
        BminusVKnownList = []
        BminusVMeasuredList = []
        BminusbList = []
        VminusvList = []
        for i in range(0, len(ObjectList)):
            objectName = ObjectList[i]
            date = dateList[i]
            stackNum = 8
            if method == 'median':
                if os.path.exists(f'Studentship/PhotutilsTutorials/{date}_{objectName}_MedianStackedApertureResults-B_{stackNum}.txt') == False:
                    imageStack(objectName, 'B', method, date, stackNum)
                    imageStack(objectName, 'V', method, date, stackNum)
                    PSFStackedPhotometry(1, 'V', objectName, method, date, stackNum)
                    annulusPhotometry(objectName, 'V', method, date, stackNum)
                    PSFStackedPhotometry(1, 'B', objectName, method, date, stackNum)
                    annulusPhotometry(objectName, 'B', method, date, stackNum)

            if method == 'mean':
                if os.path.exists(f'Studentship/PhotutilsTutorials/{date}_{objectName}_MeanStackedApertureResults-B_{stackNum}.txt') == False:
                    imageStack(objectName, 'B', method, date, stackNum)
                    imageStack(objectName, 'V', method, date, stackNum)
                    PSFStackedPhotometry(1, 'V', objectName, method, date, stackNum)
                    annulusPhotometry(objectName, 'V', method, date, stackNum)
                    PSFStackedPhotometry(1, 'B', objectName, method, date, stackNum)
                    annulusPhotometry(objectName, 'B', method, date, stackNum)
            

            file1 = f'Studentship/PhotutilsTutorials/{objectName}Stars.txt'
            if method == 'median':
                file2 = f'Studentship/PhotutilsTutorials/{date}_{objectName}_MedianStackedApertureResults-B_{stackNum}.txt'
                file3 = f'Studentship/PhotutilsTutorials/{date}_{objectName}_MedianStackedApertureResults-V_{stackNum}.txt'
            else:
                file2 = f'Studentship/PhotutilsTutorials/{date}_{objectName}_MeanStackedApertureResults-B_{stackNum}.txt'
                file3 = f'Studentship/PhotutilsTutorials/{date}_{objectName}_MeanStackedApertureResults-V_{stackNum}.txt'

            with open(file2) as f4:
                lines4 = f4.readlines()

            i = 0
            ignoreList = []
            for line4 in lines4:
                var = line4.strip().strip('\n').split(' ')
                if float(var[3]) < 0:
                    ignoreList.append(i)
                i += 1

            with open(file1, 'r') as f:
                lines = f.readlines()

            BFluxList = []
            VFluxList = []

            knownMags_B = []
            knownMags_V = []
            j = 0
            for line in lines:
                var = line.strip('\n').split(' ')
                if j not in ignoreList:
                    knownMags_B.append(float(var[2]))
                    knownMags_V.append(float(var[3]))
                j += 1

            with open(file2, 'r') as f2:
                lines2 = f2.readlines()

            for line2 in lines2:
                var2 = line2.strip('\n').split(' ')
                if float(var2[3]) > 0:
                    BFluxList.append(float(var2[3]))

            with open(file3, 'r') as f3:
                lines3 = f3.readlines()

            for line3 in lines3:
                var3 = line3.strip('\n').split(' ')
                if float(var3[3]) > 0:
                    VFluxList.append(float(var3[3]))

            print(np.shape(VFluxList))
            print(np.shape(BFluxList))
            standardMag_B = []
            standardMag_V = []
            

            knownMagB = knownMags_B[0]
            knownMagV = knownMags_V[0]
            knownFluxB = BFluxList[0]
            knownFluxV = VFluxList[0]

            i = 0
            for i in range(0, len(BFluxList)):
                targetFluxB = float(BFluxList[i])
                targetFluxV = float(VFluxList[i])
                standardMagTarget_B = knownMagB -2.5*np.log10(targetFluxB / knownFluxB)
                standardMagTarget_V = knownMagV -2.5*np.log10(targetFluxV / knownFluxV)
                standardMag_B.append(standardMagTarget_B)
                standardMag_V.append(standardMagTarget_V)
        

            standardMag_B = np.array(standardMag_B[1:]) 
            standardMag_V = np.array(standardMag_V[1:])
            knownMagsB = np.array(knownMags_B[1:])
            knownMagsV = np.array(knownMags_V[1:])

            BminusVMeasured = np.subtract(standardMag_B, standardMag_V)
            BminusVMeasuredList.append(BminusVMeasured)
            BminusVKnown = np.subtract(knownMagsB, knownMagsV)
            BminusVKnownList.append(BminusVKnown)

            Bminusb = np.subtract(knownMagsB, standardMag_B)
            BminusbList.append(Bminusb)
            Vminusv = np.subtract(knownMagsV, standardMag_V)
            VminusvList.append(Vminusv)
            
        if mode == 'single':
            paramsList1 = []
            paramsList2 = []
            paramsList3 = []    
            for j in range(0, len(BminusVKnownList)):
                newXY1 = OutlierRegress(BminusVKnownList[j], BminusVMeasuredList[j])
                paramsList1.append(linregress(newXY1[0], newXY1[1]))
                newXY2 = OutlierRegress(BminusVKnownList[j], BminusbList[j])
                paramsList2.append(linregress(newXY2[0], newXY2[1]))
                newXY3 = OutlierRegress(BminusVKnownList[j], VminusvList[j])
                paramsList3.append(linregress(newXY3[0], newXY3[1]))
            
            for jj in range(0,len(paramsList1)):
                grad = float(paramsList1[jj].slope)
                graderror = paramsList1[jj].stderr
                intercept = paramsList1[jj].intercept
                Tbv = 1/grad
                Tbv_error = graderror / (grad**2)
                x = np.array([-0.3, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
                y = grad*x + intercept
                print(f'T_bv is {Tbv} +- {Tbv_error}\n')
                plt.plot(x,y)
                plt.scatter(BminusVKnownList[jj], BminusVMeasuredList[jj], label = dateList[jj])
            plt.xlabel("Known B-V")
            plt.ylabel("Calculated B-V")
            plt.legend()
            plt.title(f'{objectName} B-V, Calculated vs Known')
            plt.show()

            for xx in range(0, len(paramsList2)):
                Tb_bv = paramsList2[xx].slope
                Tb_bv_error = paramsList2[xx].stderr
                intercept2 = paramsList2[xx].intercept
                plt.scatter(BminusVKnownList[xx], BminusbList[xx], label = dateList[xx])
                x = np.array([-0.3, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
                y = Tb_bv*x + intercept2
                plt.plot(x,y)
                print(f'Tb_bv is {Tb_bv} +- {Tb_bv_error}\n')
            plt.xlabel("Known B-V")
            plt.ylabel("Measured B-b")
            plt.legend()
            plt.title(f'{objectName} B-b Calculated vs B-V Known')
            plt.show()

            for zz in range(0, len(paramsList3)):
                Tv_bv = paramsList3[zz].slope
                Tv_bv_error = paramsList3[zz].stderr
                intercept3 = paramsList3[zz].intercept
                plt.scatter(BminusVKnownList[zz], VminusvList[zz], label = dateList[zz])
                x = np.array([-0.3, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
                y = Tv_bv*x + intercept3
                plt.plot(x,y)
                print(f'Tv_bv is {Tv_bv} +- {Tv_bv_error}\n')
            plt.xlabel("Known B-V")
            plt.ylabel("Measured V-v")
            plt.legend()
            plt.title(f'{objectName} V-v Calculated vs B-V Known')
            plt.show()

        elif mode == 'all':
            BminusVKnown = []
            BminusVMeasured = []
            Bminusb = []
            Vminusv = []
            for i in range(0, len(BminusVKnownList)):
                for j in range(0, len(BminusVKnownList[i])):
                    BminusVKnown.append(BminusVKnownList[i][j])
                    BminusVMeasured.append(BminusVMeasuredList[i][j])
                    Bminusb.append(BminusbList[i][j])
                    Vminusv.append(VminusvList[i][j])
            newXY1 = OutlierRegress(BminusVKnown, BminusVMeasured)
            params1 = linregress(newXY1[0], newXY1[1])
            newXY2 = OutlierRegress(BminusVKnown, Bminusb)
            params2 = linregress(newXY2[0], newXY2[1])
            newXY3 = OutlierRegress(BminusVKnown, Vminusv)
            params3 = linregress(newXY3[0], newXY3[1])

            
            grad = float(params1.slope)
            graderror = params1.stderr
            intercept = params1.intercept
            Tbv = 1/grad
            Tbv_error = graderror / (grad**2)
            x = np.array([-0.3, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
            y = grad*x + intercept
            print(f'T_bv is {Tbv} +- {Tbv_error}\n')
            plt.plot(x,y)
            plt.scatter(BminusVKnown, BminusVMeasured)
            plt.xlabel("Known B-V")
            plt.ylabel("Calculated B-V")
            plt.show()

            
            Tb_bv = params2.slope
            Tb_bv_error = params2.stderr
            intercept2 = params2.intercept
            plt.scatter(BminusVKnown, Bminusb)
            x = np.array([-0.3, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
            y = Tb_bv*x + intercept2
            plt.plot(x,y)
            print(f'Tb_bv is {Tb_bv} +- {Tb_bv_error}\n')
            plt.xlabel("Known B-V")
            plt.ylabel("Measured B-b")
            plt.show()

            Tv_bv = params3.slope
            Tv_bv_error = params3.stderr
            intercept3 = params3.intercept
            plt.scatter(BminusVKnown, Vminusv)
            x = np.array([-0.3, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
            y = Tv_bv*x + intercept3
            plt.plot(x,y)
            print(f'Tv_bv is {Tv_bv} +- {Tv_bv_error}\n')
            plt.xlabel("Known B-V")
            plt.ylabel("Measured V-v")
            plt.show()


#FluxCalibration('MCT0401', 'mean', '20230114')
#RepeatedFluxCalibration('JL163,JL163,JL163', '20221214,20230115,20230116', 'mean', 'single')
















