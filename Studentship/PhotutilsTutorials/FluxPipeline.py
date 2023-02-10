from projectMultiple import imageProject
from stackMultiple import imageStack
from psfForStackedImages import PSFStackedPhotometry
from FinalAperturePhotometry import annulusPhotometry
from FluxCalibration1 import FluxCalibration
import os

def pathCheck(objectName, method, date, stackNum):
    if os.path.exists(f'Studentship/PhotutilsTutorials/{date}_{objectName}_{stackNum}-Coefficients.txt'):
        os.remove(f'Studentship/PhotutilsTutorials/{date}_{objectName}_{stackNum}-Coefficients.txt')
    if os.path.exists(f'Studentship/PhotutilsTutorials/{date}_{objectName}_{stackNum}-Magnitudes.txt'):
        os.remove(f'Studentship/PhotutilsTutorials/{date}_{objectName}_{stackNum}-Magnitudes.txt')
    if os.path.exists(f'Studentship/PhotutilsTutorials/{date}_{objectName}_MeanStackedApertureResults-B_{stackNum}.txt'):
        os.remove(f'Studentship/PhotutilsTutorials/{date}_{objectName}_MeanStackedApertureResults-B_{stackNum}.txt')
        os.remove(f'Studentship/PhotutilsTutorials/{date}_{objectName}_MeanStackedApertureResults-V_{stackNum}.txt')
    if os.path.exists(f"Studentship/PhotutilsTutorials/{date}_{objectName}Stacked_B_{stackNum}_PSF.txt"):
        os.remove(f"Studentship/PhotutilsTutorials/{date}_{objectName}Stacked_B_{stackNum}_PSF.txt")
        os.remove(f"Studentship/PhotutilsTutorials/{date}_{objectName}Stacked_V_{stackNum}_PSF.txt")
    if os.path.exists(f"Studentship/PhotutilsTutorials/{date}_{objectName}Stacked_B_{stackNum}_mean_PSF.txt"):
        os.remove(f"Studentship/PhotutilsTutorials/{date}_{objectName}Stacked_B_{stackNum}_mean_PSF.txt")
        os.remove(f"Studentship/PhotutilsTutorials/{date}_{objectName}Stacked_V_{stackNum}_mean_PSF.txt")
    if os.path.exists(f"Studentship/PhotutilsTutorials/{date}_{objectName}_MedianStackedApertureResults-B_{stackNum}.txt"):
        os.remove(f"Studentship/PhotutilsTutorials/{date}_{objectName}_MedianStackedApertureResults-B_{stackNum}.txt")
        os.remove(f"Studentship/PhotutilsTutorials/{date}_{objectName}_MedianStackedApertureResults-V_{stackNum}.txt")


def automatedAnalysis(objectName, File1, File2, method, date, stackNum, mode, bkgsub, overwrite:bool, increment):
    """
    The idea is that this function will do both median and mean stacking for a single object, and will 
    create incrementally stacked images according to the specified increment you specify.\n
    \n
    objectName - string that contains object name, e.g. 'JL163 '\n
    File1 - File path of base image upon which all others are reprojected, B Filter \n
    File2 - same but V filter \n
    method - 'median', 'mean' or 'both' \n
    date - string of all dates separated by commas, in yyyymmdd format. e.g. '20221214,20230112,20230114...' \n
    stackNum - 8 times the amount of dates, usually, otherwise it will be the maximum stack size created \n
    mode - always use multiple \n
    bkgsub - if doing meaningful data analysis, must be True \n
    overwrite - True if you want to overwrite any existing files. False otherwise. \n

    """
    startVal = stackNum
    originalIncrement = increment
    while startVal > 0:
        if method == 'both':
            FluxCalibrationAnalysis(objectName, File1, File2, 'mean', date, startVal, mode, bkgsub, overwrite)
            FluxCalibrationAnalysis(objectName, File1, File2, 'median', date, startVal, mode, bkgsub, overwrite)
        elif method == 'mean':
            FluxCalibrationAnalysis(objectName, File1, File2, 'mean', date, startVal, mode, bkgsub, overwrite)
        elif method == 'median':
            FluxCalibrationAnalysis(objectName, File1, File2, 'median', date, startVal, mode, bkgsub, overwrite)
        startVal -= originalIncrement
        

def FluxCalibrationAnalysis(objectName, File1, File2, method, date, stackNum, mode, bkgsub, overwrite:bool):

    """
    Should almost always be run with stackNum equal to the maximum amount of images for that
    given dataset. Then you can run the calibration again later on the smaller levels.

    objectName - string that contains object name, e.g. 'JL163 '\n
    File1 - File path of base image upon which all others are reprojected, B Filter \n
    File2 - same but V filter \n
    method - 'median', 'mean' or 'both' \n
    date - string of all dates separated by commas, in yyyymmdd format. e.g. '20221214,20230112,20230114...' \n
    stackNum - 8 times the amount of dates, usually, otherwise it will be the maximum stack size created \n
    mode - always use multiple \n
    bkgsub - if doing meaningful data analysis, must be True \n
    overwrite - True if you want to overwrite any existing files. False otherwise. \n
    """

    ## ADD A CHECKING MECHANISM TO CHECK THAT CERTAIN REPROJECTED IMAGES ALREADY EXIST
    ## AS TO NOT REPEAT THE REPROJECTION PROCESS. CAN ALSO DO FOR STACKING
    path1 = f'Studentship/ReprojectedImages/{date}_{objectName}Reprojected-{stackNum}-B.fits'
    path2 = f'Studentship/ReprojectedImages/{date}_{objectName}Reprojected-{stackNum}-V.fits'
    path3 = f'Studentship/StackedImages/{date}_{objectName}Stacked_mean_B_{stackNum}.fits'
    path4 = f'Studentship/StackedImages/{date}_{objectName}Stacked_mean_V_{stackNum}.fits'
    path5 = f'Studentship/StackedImages/{date}_{objectName}Stacked_{stackNum}_B.fits'
    path6 = f'Studentship/StackedImages/{date}_{objectName}Stacked_{stackNum}_V.fits'
    pathCheck(objectName, method, date, stackNum)

    if overwrite == True:
        if method == 'mean':
            if os.path.exists(f'Studentship/StackedImages/{date}_{objectName}Stacked_mean_V_{stackNum}.fits'):
                os.remove(f'Studentship/StackedImages/{date}_{objectName}Stacked_mean_B_{stackNum}.fits')
                os.remove(f'Studentship/StackedImages/{date}_{objectName}Stacked_mean_V_{stackNum}.fits')
        elif method == 'median':
            if os.path.exists(f'Studentship/StackedImages/{date}_{objectName}Stacked_{stackNum}_B.fits'):
                os.remove(f'Studentship/StackedImages/{date}_{objectName}Stacked_{stackNum}_B.fits')
                os.remove(f'Studentship/StackedImages/{date}_{objectName}Stacked_{stackNum}_V.fits')

    if (os.path.exists(path1) == False) and (os.path.exists(path2) == False):
        imageProject(objectName, 'V', File1, date, mode)
        imageProject(objectName, 'B', File2, date, mode)
    if (overwrite == True) or (((os.path.exists(path3) == False) and (os.path.exists(path4) == False)) or ((os.path.exists(path5) == False) and (os.path.exists(path6) == False))):
        imageStack(objectName, 'V', method, date, bkgsub, stackNum)
        imageStack(objectName, 'B', method, date, bkgsub, stackNum)
    PSFStackedPhotometry(1, 'V', objectName, method, date, stackNum)
    annulusPhotometry(objectName, 'V', method, date, stackNum)
    PSFStackedPhotometry(1, 'B', objectName, method, date, stackNum)
    annulusPhotometry(objectName, 'B', method, date, stackNum)

    FluxCalibration(objectName, method, date, stackNum)
    

#FluxCalibrationAnalysis('JL163', 'Studentship/SPIRIT6-Observations/GoodDarks/Multiple-Exposures/20221214/JL-163/JL-163-S001-V-R001-V.fts', 'Studentship/SPIRIT6-Observations/GoodDarks/Multiple-Exposures/20221214/JL-163/JL-163-S001-B-R001-B.fts', 'mean', '20221214,20230115,20230116', 10, 'multiple', True, True) #16 is done, 14 then 12 are next
automatedAnalysis('JL163', 'Studentship/SPIRIT6-Observations/GoodDarks/Multiple-Exposures/20221214/JL-163/JL-163-S001-V-R001-V.fts', 'Studentship/SPIRIT6-Observations/GoodDarks/Multiple-Exposures/20221214/JL-163/JL-163-S001-B-R001-B.fts', 'both', '20221214,20230115,20230116', 2, 'multiple', True, True, 2)