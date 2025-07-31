# Title: Free shifting algorithms. 
#
# Author: Farhad Abdollahi (farhad.abdollahi.ctr@dot.gov)
# Date: 01/28/2025
# ======================================================================================================================

# Importing the required libraries.
import os
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar, curve_fit, minimize
from scipy.interpolate import Akima1DInterpolator


def FreeShift_GordonShaw1994(data, Tref, Tglassy, DensityShift):
    """
    This function performs the Gordon-Shaw free shifting algorithm. It will freely shift the frequency sweep test 
    isotherms to produce a smooth shifted curve. 
    It is noted that this function and its sub-functions are based on the FORTRAN code, originally developed by 
    Gordon & Shaw (1994) with a little modification during the translation to Python. Reference:
        Gordon, Shaw (1994), Computer programs for rheologists, Hansen.

    :param data: A dataframe of all isotherms.
    :param Tref: Reference temperature (C).
    :param Tglassy: Glassy temperature as an input (C).
    :param DensityShift: A dictionary to whether perform the density shift or not. If so, the required parameters 
    included. 
    """
    # Use only the non-outlier datapoints.
    Tempdf = data[data['IsOutlier'] == 0]
    # Find the unique temperature isotherms. 
    UniqueTemp = np.sort(Tempdf['Temperature'].unique())       # all different temperatures. 
    n_temp = len(UniqueTemp)
    # Extract the datapoints, for easier use.
    Temp    = Tempdf['Temperature'].to_numpy()
    Freq    = Tempdf['Frequency'].to_numpy()
    Gstar   = Tempdf['|G*|'].to_numpy()
    Phase   = Tempdf['PhaseAngle'].to_numpy()
    Storage = Gstar * np.cos(np.deg2rad(Phase))
    Loss    = Gstar * np.sin(np.deg2rad(Phase))
    # Prepare some variables for analysis (Keep the FORTRAN structure). 
    Shift_cell = [[], [], []]
    Block1 = {'AORIG': 1, 'X': [], 'Y': [], 'I': 0, 'ND': [], 'NF': np.zeros(n_temp).astype(int), 'NP': 0}
    Block2 = {'AT': np.zeros(n_temp), 'ATTOL': [], 'T': UniqueTemp, 'TG': 0, 'TR': 0, 'NT': n_temp}
    Block2['ATTOL'] = 4                 # Acceptable number of decades for wild shift
    Block2['TR'] = Tref                 # Reference temperature as defined by user. 
    Block2['TG'] = Tglassy              # Glassy temperature as input. 
    IFLAG = 2  # Always true in the provided code (Availability of loss modulus data)
    ITYPE = 1  # Frequency scale (2 would be for time scale, which is not implemented here)
    # ------------------------------------------------------------------------------------------------------------------
    # OPTIONAL: Density-temerature correction:
    #   TFAC parameter includes the density shift which will remain zero when density shifting is not active.
    TFAC = np.zeros(n_temp)
    if DensityShift['Active']:
        for i in range(n_temp):
            TFAC[i] = np.log10(((Block2['TR'] + 273.15) / (Block2['T'][i] + 273.15)) / 
                                RHOR(Block2['T'][i] + 273.15, Block2['TG'] + 273.15, Block2['TR'] + 273.15, 
                                     DensityShift['a'][0], DensityShift['a'][1]))
    # Fill some new variables with frequency, storage, and loss moduli (Keep FORTRAN NAMES).
    F    = [[] for _ in range(Block2['NT'])]
    GP   = [[] for _ in range(Block2['NT'])]
    GPP  = [[] for _ in range(Block2['NT'])]
    for i in range(n_temp):
        Index = np.where(Temp == UniqueTemp[i])[0]
        F[i]    = np.log10(Freq[Index])
        GP[i]   = np.log10(Storage[Index] + TFAC[i])
        GPP[i]  = np.log10(Loss[Index] + TFAC[i])
        Block1['NF'][i] = len(F[i])
    # ------------------------------------------------------------------------------------------------------------------
    # Perform shifting
    IMC = 1  # Flag for Monte Carlo simulation (2 means active, not used here)
    Shift_cell[0], Block1, Block2 = SHIFT(ITYPE, IMC, F, GP,  Block1, Block2)
    Shift_cell[1], Block1, Block2 = SHIFT(ITYPE, IMC, F, GPP, Block1, Block2)
    Shift_cell[2] = (np.array(Shift_cell[0]) + np.array(Shift_cell[1])) / 2
    AT = Shift_cell[2]      # Use the average shift factors onwards. 
    # ------------------------------------------------------------------------------------------------------------------
    # # Fit a smoothing curve to G' and G" to balance the shift factors for both of them. 
    # #                       log G' = a.log(x)^10 + b.log(x)^9+...+ j.log(x) + k
    # # Calculating the reduced frequencies based on the loss and storage moduli.
    # reduced_s, reduced_l = [], []   # To store the reduced frequencies based on the storage and loss moduli. 
    # Storage, Loss = [], []          # To store the storage and loss moduli values. 
    # ATs = Shift_cell[0]             # Shift factors based on the storage modulus. 
    # ATl = Shift_cell[1]             # Shift factors based on the loss modulous. 
    # for i in range(n_temp):
    #     reduced_s.append(F[i] * (10 ** ATs[i]))
    #     reduced_l.append(F[i] * (10 ** ATl[i]))
    #     Storage.append(GP[i])
    #     Loss.append(GPP[i])
    # # Balancing the reduced frequency and moduluses.
    # # Fitting section.
    # NumData_s = len(reduced_s)
    # NumData_l = len(reduced_l)
    # LogReduced_s = np.log10(np.array(reduced_s))
    # LogReduced_l = np.log10(np.array(reduced_l))
    # LogStorage   = np.log10(np.array(Storage))
    # LogLoss      = np.log10(np.array(Loss))
    # Xs = np.column_stack([np.ones(NumData_s)] + [LogReduced_s ** p for p in range(10, 0, -1)])
    # Xl = np.column_stack([np.ones(NumData_l)] + [LogReduced_l ** p for p in range(10, 0, -1)])
    # Beta_s = np.linalg.inv(Xs.T @ Xs) @ (Xs.T @ LogStorage)
    # Beta_l = np.linalg.inv(Xl.T @ Xl) @ (Xl.T @ LogLoss)
    # # Define the functions. 
    # Fit_Func_s = lambda x: 10 ** (Beta_s[0]  + 
    #                               Beta_s[1]  * (np.log10(x) ** 10) + 
    #                               Beta_s[2]  * (np.log10(x) ** 9) + 
    #                               Beta_s[3]  * (np.log10(x) ** 8) + 
    #                               Beta_s[4]  * (np.log10(x) ** 7) + 
    #                               Beta_s[5]  * (np.log10(x) ** 6) + 
    #                               Beta_s[6]  * (np.log10(x) ** 5) + 
    #                               Beta_s[7]  * (np.log10(x) ** 4) + 
    #                               Beta_s[8]  * (np.log10(x) ** 3) + 
    #                               Beta_s[9]  * (np.log10(x) ** 2) + 
    #                               Beta_s[10] * (np.log10(x) ** 1))
    # Fit_Func_l = lambda x: 10 ** (Beta_l[0]  + 
    #                               Beta_l[1]  * (np.log10(x) ** 10) + 
    #                               Beta_l[2]  * (np.log10(x) ** 9) + 
    #                               Beta_l[3]  * (np.log10(x) ** 8) + 
    #                               Beta_l[4]  * (np.log10(x) ** 7) + 
    #                               Beta_l[5]  * (np.log10(x) ** 6) + 
    #                               Beta_l[6]  * (np.log10(x) ** 5) + 
    #                               Beta_l[7]  * (np.log10(x) ** 4) + 
    #                               Beta_l[8]  * (np.log10(x) ** 3) + 
    #                               Beta_l[9]  * (np.log10(x) ** 2) + 
    #                               Beta_l[10] * (np.log10(x) ** 1))
    # # Calculate the average reduced frequency. 
    # RedFreqNew = 10 ** ((LogReduced_s + LogReduced_l) / 2)
    # StorageNew = Fit_Func_s(RedFreqNew)
    # LossNew    = Fit_Func_l(RedFreqNew)
    # GstarNew   = np.sqrt(StorageNew ** 2 + LossNew ** 2)
    # PhaseNew   = np.rad2deg(np.arctan(LossNew / StorageNew))
    # ------------------------------------------------------------------------------------------------------------------
    # Calculate the reduced frequency and add it to the "data" dataframe.
    ShiftMapping = dict(zip(UniqueTemp, AT))
    data['ShiftFactor']   = data['Temperature'].map(ShiftMapping)
    data['Red_Frequency'] = data['Frequency'] * (10 ** data['ShiftFactor'])
    # Return the calculated shift factors, and the updated "data" dictionary. 
    return ShiftMapping, data
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


def SHIFT(ITYPE, IMC, F, GP, Block1, Block2):
    """
    Perform the superposition procedure of Gordon and Shaw method.
    """
    # Define a variable ERRAT
    ERRAT = np.zeros(51)
    # Loop over each unique temperature. 
    for i in range(Block2['NT'] - 1, 0, -1):
        Block1['I'] = i
        Block1['NP'] = Block1['NF'][i] + Block1['NF'][i - 1]
        Block1['X'] = np.concatenate((F[i], F[i - 1]))
        Block1['Y'] = np.concatenate((GP[i], GP[i - 1]))
        Coeff = np.polyfit(Block1['X'], Block1['Y'], 1)
        Block2['AT'][i] = -(np.mean(GP[i]) - np.mean(GP[i - 1])) / Coeff[0]
        # Estimate shift using the Kaelble modification WLF equation.
        ATWLF = (-1) ** (ITYPE + 1) * \
            (17.4 * (Block2['T'][i] - Block2['TG']) / (51.6 + abs(Block2['T'][i] - Block2['TG'])) - \
             17.4 * (Block2['T'][i - 1] - Block2['TG']) / (51.6 + abs(Block2['T'][i - 1] - Block2['TG'])))
        # Check for incorrect sign of shift due to widely spaced data.
        if Block2['AT'][i] * ATWLF < 0:
            Block2['AT'][i] = ATWLF
        Block1['ND'] = int(min([(len(Block1['X']) ** 0.333) + 2, 
                                max([abs(F[i - 1][-1] - F[i - 1][0]) + 2, 
                                     abs(F[i][-1] - F[i][0]) + 2, 
                                     abs(ATWLF),]),]))
        # Set up CALL to IMSL DUVMIF (using 'fminbnd' in MATLAB, or 'minimize_scalar' in Python).
        AGUESS = Block2['AT'][i]
        BOUND  = Block2['ATTOL']
        ATI    = minimize_scalar(lambda at: FAT1(at, Block1), 
                                 bounds=(AGUESS - BOUND, AGUESS + BOUND), method="bounded").x
        FX     = FAT1(ATI, Block1)
        # Estimate error in Shift factor
        ScaleFactor = 1.0           # Scaling factor for computing finite differences.
        _, _, COV = Error(lambda at: FAT1(at, Block1), ATI, ScaleFactor, FX)
        
        # Check for wild shift.
        if np.abs(ATI - Block2['AT'][i]) < Block2['ATTOL']:
            Block2['AT'][i] = ATI
            if IMC == 1:
                # Estimate error in the shift factor ERRAT.
                ERRAT[i] = T95(Block1['NP'] - 1) * np.sqrt(FX / (Block1['NP'] - 1)) * COV
        else:
            # Assume error scales with shift.
            if IMC == 1:
                ERRAT[i] = ERRAT[i+1] * np.exp(abs(ATI)) * np.exp(np.abs(Block2['AT'][i]))
    # ------------------------------------------------------------------------------------------------------------------
    # Switch shift to lower temperature of each pair.
    for i in range(Block2['NT'] - 1):
        if IMC == 1:
            ERRAT[i] = ERRAT[i+1]
        Block2['AT'][i] = Block2['AT'][i+1]
    Block2['AT'][Block2['NT'] - 1] = 0
    if IMC == 1:
        ERRAT[Block2['NT']] = ERRAT[Block2['NT'] + 1]
    for i in range(Block2['NT'] - 2, -1, -1):
        if IMC == 1:
            ERRAT[i] = ERRAT[i] + ERRAT[i + 1]
        Block2['AT'][i] = Block2['AT'][i] + Block2['AT'][i + 1]
    # ------------------------------------------------------------------------------------------------------------------
    # Interpolate to get shift at reference temperature TR.
    AT_Akima_Func = Akima1DInterpolator(Block2['T'], Block2['AT'])
    ATR = AT_Akima_Func(Block2['TR'])
    if IMC ==1:
        ERRAT_Akima_Func = Akima1DInterpolator(Block2['T'], ERRAT[:len(Block2['T'])])
        ERRATR = ERRAT_Akima_Func(Block2['TR'])
    # ------------------------------------------------------------------------------------------------------------------
    for i in range(Block2['NT']):
        if IMC == 1:
            ERRAT[i] = np.abs(ERRAT[i] - ERRATR)
        Block2['AT'][i] = Block2['AT'][i] - ATR
    AT = Block2['AT']
    # ------------------------------------------------------------------------------------------------------------------
    # Return the results. 
    return AT, Block1, Block2
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


def FAT(N, AT, Block1):
    """
    Evaluate the function (Shift factor) to be minimized.

    :param N: Degree of polynomial.
    :param AT: _description_
    :return: _description_
    """
    # Defining the parameters.
    XT = np.zeros(Block1['NP'])
    W = np.zeros(Block1['NP'])

    # Check for the sign of AT to be correct.
    if (Block1['AORIG'] * AT) < 0:
        Ysd = Block1['Y'].std()
        F = (Ysd ** 2) * (len(Block1['Y']) - N)
        return F
    
    # Generating the XT array from initially shifted frequencies.
    for j in range(int(Block1['NF'][Block1['I']])):
        XT[j] = Block1['X'][j]
    for j in range(int(Block1['NF'][Block1['I']]), Block1['NP']):
        XT[j] = Block1['X'][j] + AT

    # Get residuals R(NP) from polynomial fit and weighted factors, w(NP).
    _, _, R, _ = POLYN(XT, Block1['Y'], Block1['ND'])
    XAVE  = XT.mean()
    XWIDE = (XT.max() - XT.min()) / 2
    # estimate the Lorentzian weighted factors for residuals.
    for j in range(Block1['NP']):
        W[j] = 1 / (1 + ((XT[j] - XAVE) / XWIDE) ** 2)

    # Calculating the result.
    F = 0
    for j in range(Block1['NP']):
        F = F + W[j] * (R[j] ** 2)
    
    # Return the result.
    return F
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


def FAT1(AT, Block1):
    """
    This function compute the function value (Shift factor) to be minimized via function FAT. 

    :param AT: _description_
    :return: _description_
    """
    return FAT(1, AT, Block1)
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


def T95(N):
    """
    t-statistics for 95% confidence interval with N degrees of freedom.
    """
    return 1.96 + 2.3539 / N + 3.1902 / (N ** 2) + 0.7004 / (N ** 3) + 4.5015 / (N ** 4)
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


def RHOR(T, TG, TR, AG, AR):
    """
    Calculate the density temperature correction.
    """
    if T >= TG and TR >= TG:
        return np.exp(-AR * (T - TR))
    elif T < TG and TR < TG:
        return np.exp(-AG * (T - TR))
    elif T >= TG and TR < TG:
        return np.exp(-AR * (T - TG)) * np.exp(-AG * (TG - TR))
    elif T < TG and TR >= TG:
        return np.exp(-AR * (TG - TR)) * np.exp(-AG * (T - TR))
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


def POLYN(X, Y, N):
    """
    Fit a polynomial of (N-1)th order to the X and Y data.

    Parameters:
        X (array-like): The input data for the independent variable.
        Y (array-like): The input data for the dependent variable.
        N (int): The degree of the polynomial + 1.

    Returns:
        B (numpy.ndarray): Coefficients of the polynomial.
        ERR (numpy.ndarray): Error matrix (equivalent to covariance matrix of the fit).
        RES (numpy.ndarray): Residuals.
        SEE (float): Standard error of estimation.
    """
    n_data = len(X)
    # Fit polynomial of degree N-1
    B, residuals, _, _, _ = np.polyfit(X, Y, N-1, full=True)
    # Residuals
    RES = Y - np.polyval(B, X)
    # Sum of squared residuals
    sum_squared_residuals = np.sum(RES ** 2)
    # Standard error of estimation
    SEE = np.sqrt(sum_squared_residuals / (n_data - N))
    # Error matrix (covariance matrix approximation)
    # Placeholder for ERR (not directly available in numpy)
    # You can use additional libraries like scipy.stats for detailed error estimates
    ERR = np.cov(X, Y) if len(residuals) > 0 else None
    # Return the results. 
    return B, ERR, RES, SEE
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


def Error(fcn, ati, ascale, fx):
    """
    This function calculates the diagonal elements of the covariance matrix (COV) estimated from the curvature 
    matrix (HI). It approximates the Hessian matrix (H) using forward differences at the function value fx.

    Parameters:
    - fcn: Function to evaluate.
    - ati: Input parameter at which to evaluate the function.
    - ascale: Scaling factor for computing finite differences.
    - fx: Function value at ati (not used in this code).

    Returns:
    - H: Approximated Hessian matrix.
    - HI: Inverse of the curvature matrix.
    - COV: Covariance vector (equal to HI in this case).
    """
    EPSFCN = 0.22e-15
    # Compute step size for finite differences
    h = (EPSFCN ** (1 / 3)) * max(abs(ati), 1 / ascale) * np.sign(ati)
    # Calculate the Hessian matrix using finite differences
    H = (-1 * fcn(ati + 2 * h) 
         + 16 * fcn(ati + h) 
         - 30 * fcn(ati) 
         + 16 * fcn(ati - h) 
         - 1 * fcn(ati - 2 * h)) / (12 * h ** 2)
    # Calculate the curvature matrix
    H = H / 2
    # Inverse of the curvature matrix
    HI = 1 / H
    # Covariance vector
    COV = HI
    # Return the results. 
    return H, HI, COV
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


def ShiftModel_WLF(Temp, Tref, C1, C2): 
    """
    This function calculates the shift factor using the Williams-Landel-Ferry (WLF) model. 

    :param Temp: An array of the "Temperatures" (°C).
    :param Tref: Reference temperature (°C).
    :param C1: Model coefficient.
    :param C2: Model coefficient.
    """
    return -C1 * (Temp - Tref) / (C2 + Temp - Tref)
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


def Fit_ShiftModel_WLF(Temp, aT, Tref, InitialGuess=[10.0, 55.0]):
    """
    This function fits the Williams-Landel-Ferry (WLF) model to the calculated shift factors at different temperatures. 
    It only used for fitting a WLF model to the results of Free Shifting, for comparison purposes. 

    :param Temp: A list of the isotherm temperatures (°C).
    :param aT: A list of the shift factors.
    :param Tref: The reference temperature (°C)
    """
    model = lambda T, C1, C2: ShiftModel_WLF(T, Tref, C1, C2)                       # Define the WLF model. 
    FitCoeff, _ = curve_fit(model, np.array(Temp), np.array(aT), p0=InitialGuess)   # Fit the model. 
    C1, C2 = FitCoeff
    # Return the results. 
    return C1, C2
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


def Fit_MC_ChristensenAnderson(data, FitType, InitialGuess=[]):
    """
    This function fits the Christensen-Anderson (CA) model to the already shifted datapoints. It can fit to |G*|, G', 
    G", or δ. 

    :param data: A dataframe of the datapoints. 
    :param FitType: An integer for fit type, 0: to |G*|, 1: to δ, 2: to G', and 3: to G".
    :param InitialGuess: A list of initial guesses for the model parameters. 
    :return: A dictionary of the fitted coefficients, as well as the error metrics. 
    """
    # Use only the non-outlier datapoints.
    Tempdf = data[data['IsOutlier'] == 0]
    # Extract the data.
    RedFreq = Tempdf['Red_Frequency'].to_numpy()
    Gstar   = Tempdf['|G*|'].to_numpy()
    Phase   = Tempdf['PhaseAngle'].to_numpy()
    Storage = Tempdf['StorageModulus'].to_numpy()
    Loss    = Tempdf['LossModulus'].to_numpy()
    Freq    = np.logspace(np.log10(RedFreq.min()) - 0.3, 
                          np.log10(RedFreq.max()) + 0.3, num=250)  # reduced frequencies for plotting
    NumData = len(RedFreq)
    # Define the |G*| and δ models. 
    modelGstar = lambda w, Gg, wc, R: Gg - (R / np.log10(2)) * np.log10(1 + ((10 ** wc) / w) ** (np.log10(2) / R))
    modelPhase = lambda w, wc, R: 90 / (1 + (w / (10 ** wc)) ** (np.log10(2) / R))
    # ----------------------------------------------
    # Define the weights to the datapoints. 
    Sigma = np.ones(len(RedFreq))               # Equal weights. 
    # Min, Max = np.log10(RedFreq.min()), np.log10(RedFreq.max())
    # Mid = 0.5 * (Min + Max)
    # Sigma = 1 / (1 + 1.0 * (1 - np.exp(-((np.log10(RedFreq) - Mid) ** 2 / (Max - Min)))))
    # ----------------------------------------------
    # Check the fit type. 
    if FitType == 0:        # Fit to the complex modulus. 
        # Use default initial Guess, if not provided. 
        if InitialGuess == []: 
            InitialGuess = [8.5, 2.0, 2.5]
        # Using the "Curve fit" function, which uses unbound "Levenberg-Marquardt" algorithm and non-linear least 
        #   square fitting loss function.
        popt, pcov = curve_fit(modelGstar, RedFreq, np.log10(Gstar), p0=InitialGuess, absolute_sigma=True, sigma=Sigma)
        # --------------------------------------------------------------------------------------------------------------
        # # For testing, using more general fitting, with Mean Absolute Error loss function. 
        # result = minimize(LossFuncMAE, InitialGuess, args=(RedFreq, np.log10(Gstar)), method="Nelder-Mead")
        # popt   = result.x
        # --------------------------------------------------------------------------------------------------------------
        Coeff = [popt[0], popt[1], popt[2]]
    # ------------------------------------------------------------------------------------------------------------------
    elif FitType == 1:      # Fit to the Phase angle. 
        # Use default initial Guess, if not provided. 
        if InitialGuess == []: 
            InitialGuess = [8.5, 2.0, 2.5]
        # First fit the Phase angle model with wc and R only. Using the "Curve fit" function, which uses unbound 
        #   "Levenberg-Marquardt" algorithm and non-linear least square fitting loss function.
        Ppopt, _ = curve_fit(modelPhase, RedFreq, Phase, p0=InitialGuess[1:], absolute_sigma=True, sigma=Sigma)
        # Now, defining a new model to adjust the Gg value. 
        modelGstar2 = lambda w, Gg: Gg - (Ppopt[1] / np.log10(2)) * \
            np.log10(1 + ((10 ** Ppopt[0]) / w) ** (np.log10(2) / Ppopt[1]))
        Gpopt, _ = curve_fit(modelGstar2, RedFreq, np.log10(Gstar), 
                             p0=InitialGuess[0], absolute_sigma=True, sigma=Sigma)
        Coeff = [Gpopt[0], Ppopt[0], Ppopt[1]]
    # ------------------------------------------------------------------------------------------------------------------
    elif FitType == 2:      # Fit to the Storage modulus data (G').
        # Define the storage modulus model. 
        def modelGP(w, Gg, wc, R):
            G = modelGstar(w, Gg, wc, R)
            P = modelPhase(w, wc, R)
            return np.log10((10 ** G) * np.cos(np.deg2rad(P))) 
        # Use default initial Guess, if not provided. 
        if InitialGuess == []: 
            InitialGuess = [8.5, 2.0, 2.5]
        # ----------------------------------------------
        # Using the "Curve fit" function, which uses unbound "Levenberg-Marquardt" algorithm and non-linear least 
        #   square fitting loss function.
        popt, _ = curve_fit(modelGP, RedFreq, np.log10(Storage), p0=InitialGuess, absolute_sigma=True, sigma=Sigma)
        # --------------------------------------------------------------------------------------------------------------
        # # For testing, using more general fitting, with Mean Absolute Error loss function. 
        # result = minimize(LossFuncMAE, InitialGuess, args=(RedFreq, np.log10(Gstar)), method="Nelder-Mead")
        # popt   = result.x
        # --------------------------------------------------------------------------------------------------------------
        Coeff = [popt[0], popt[1], popt[2]]
    # ------------------------------------------------------------------------------------------------------------------
    elif FitType == 3:      # Fit to the Loss modulus data (G").
        # Define the storage modulus model. 
        def modelGPP(w, Gg, wc, R):
            G = modelGstar(w, Gg, wc, R)
            P = modelPhase(w, wc, R)
            return np.log10((10 ** G) * np.sin(np.deg2rad(P))) 
        # Use default initial Guess, if not provided. 
        if InitialGuess == []: 
            InitialGuess = [8.5, 2.0, 2.5]
        # ----------------------------------------------
        # Using the "Curve fit" function, which uses unbound "Levenberg-Marquardt" algorithm and non-linear least 
        #   square fitting loss function.
        popt, _ = curve_fit(modelGPP, RedFreq, np.log10(Loss), p0=InitialGuess, absolute_sigma=True, sigma=Sigma)
        # --------------------------------------------------------------------------------------------------------------
        # # For testing, using more general fitting, with Mean Absolute Error loss function. 
        # result = minimize(LossFuncMAE, InitialGuess, args=(RedFreq, np.log10(Gstar)), method="Nelder-Mead")
        # popt   = result.x
        # --------------------------------------------------------------------------------------------------------------
        Coeff = [popt[0], popt[1], popt[2]]
    # ------------------------------------------------------------------------------------------------------------------
    else: 
        raise Exception(f'Fit type is not recognized!! Please check the code manually!')
    # ------------------------------------------------------------------------------------------------------------------
    # Prepare outputs. 
    # Calculate the MAE loss function value. 
    MLAE_Gstar = np.abs(np.log10(Gstar) - modelGstar(RedFreq, Coeff[0], Coeff[1], Coeff[2])).sum() / NumData
    MLAE_Phase = np.abs(Phase - modelPhase(RedFreq, Coeff[1], Coeff[2])).sum() / NumData
    MLAE_GP    = np.abs(np.log10(Storage) - 
                        np.log10(10 ** modelGstar(RedFreq, Coeff[0], Coeff[1], Coeff[2]) * 
                                    np.cos(np.deg2rad(modelPhase(RedFreq, Coeff[1], Coeff[2]))))).sum() / NumData
    MLAE_GPP   = np.abs(np.log10(Loss) - 
                        np.log10(10 ** modelGstar(RedFreq, Coeff[0], Coeff[1], Coeff[2]) * 
                                    np.sin(np.deg2rad(modelPhase(RedFreq, Coeff[1], Coeff[2]))))).sum() / NumData
    # Extract the results and the model predictions. 
    Gg, wc, R = 10 ** Coeff[0], 10 ** Coeff[1], Coeff[2]
    PlotGstar = 10 ** modelGstar(Freq, Coeff[0], Coeff[1], Coeff[2])
    PlotPhase = modelPhase(Freq, Coeff[1], Coeff[2])
    PlotGP    = PlotGstar * np.cos(np.deg2rad(PlotPhase))
    PlotGPP   = PlotGstar * np.sin(np.deg2rad(PlotPhase))
    # Store the results. 
    Res = {'MC_Coeff': [Gg, wc, R],
            'MC_Coeff_Labels': ['Glassy Modulus (psi)', 'Crossover Frequency (rad/s)', 'Rheology Index'],
            'Plot_Freq': Freq, 'Plot_Gstar': PlotGstar, 'Plot_Phase': PlotPhase, 
            'Plot_Storage': PlotGP, 'Plot_Loss': PlotGPP,
            'MC_Index': [R, wc, 10 ** modelGstar(wc, Coeff[0], Coeff[1], Coeff[2])], 
            'MC_Index_Labels': ['Rheology Index', 'Crossover Frequency (rad/s)', 'Crossover Modulus (psi)'], 
            'MC_Error': [MLAE_Gstar, MLAE_Phase, MLAE_GP, MLAE_GPP], 
            'MC_Error_Labels': ['MLAE (|G*|)', 'MLAE (δ)', "MLAE (G')", 'MLAE (G")']}
    # Return the results. 
    return Res
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


def Fit_MC_ChristensenAndersonMarasteanu(data, FitType, InitialGuess=[]):
    """
    This function fits the Christensen-Anderson-Marasteanu (CAM) model to the already shifted datapoints. It can fit to 
    |G*|, G', G", or δ. 

    :param data: A dataframe of the datapoints. 
    :param FitType: An integer for fit type, 0: to |G*|, 1: to δ, 2: to G', and 3: to G".
    :param InitialGuess: A list of initial guesses for the model parameters. 
    :return: A dictionary of the fitted coefficients, as well as the error metrics. 
    """
    # Use only the non-outlier datapoints.
    Tempdf = data[data['IsOutlier'] == 0]
    # Extract the data.
    RedFreq = Tempdf['Red_Frequency'].to_numpy()
    Gstar   = Tempdf['|G*|'].to_numpy()
    Phase   = Tempdf['PhaseAngle'].to_numpy()
    Storage = Tempdf['StorageModulus'].to_numpy()
    Loss    = Tempdf['LossModulus'].to_numpy()
    Freq    = np.logspace(np.log10(RedFreq.min()) - 0.3, 
                          np.log10(RedFreq.max()) + 0.3, num=250)  # reduced frequencies for plotting
    NumData = len(RedFreq)
    # Define the |G*| and δ models. 
    modelGstar = lambda w, Gg, wc, m, k: Gg - (m / k) * np.log10(1 + ((10 ** wc) / w) ** k)
    modelPhase = lambda w, wc, m, k: 90 * m / (1 + (w / (10 ** wc)) ** k)
    # ----------------------------------------------
    # Define the weights to the datapoints. 
    Sigma = np.ones(len(RedFreq))               # Equal weights. 
    # Min, Max = np.log10(RedFreq.min()), np.log10(RedFreq.max())
    # Mid = 0.5 * (Min + Max)
    # Sigma = 1 / (1 + 1.0 * (1 - np.exp(-((np.log10(RedFreq) - Mid) ** 2 / (Max - Min)))))
    # ----------------------------------------------
    # Check the fit type. 
    if FitType == 0:        # Fit to the complex modulus. 
        # Use default initial Guess, if not provided. 
        if InitialGuess == []: 
            InitialGuess = [8.5, 2.0, 1.0, 0.4]
        # Using the "Curve fit" function, which uses unbound "Levenberg-Marquardt" algorithm and non-linear least 
        #   square fitting loss function.
        popt, pcov = curve_fit(modelGstar, RedFreq, np.log10(Gstar), p0=InitialGuess, absolute_sigma=True, sigma=Sigma)
        Coeff = [popt[0], popt[1], popt[2], popt[3]]
    # ------------------------------------------------------------------------------------------------------------------
    elif FitType == 1:      # Fit to the Storage modulus data (G').
        # Define the storage modulus model. 
        def modelGP(w, Gg, wc, m, k):
            G = modelGstar(w, Gg, wc, m, k)
            P = modelPhase(w, wc, m, k)
            return np.log10((10 ** G) * np.cos(np.deg2rad(P))) 
        # Use default initial Guess, if not provided. 
        if InitialGuess == []: 
            InitialGuess = [8.5, 2.0, 1.0, 0.4]
        # Using the "Curve fit" function, which uses unbound "Levenberg-Marquardt" algorithm and non-linear least 
        #   square fitting loss function.
        popt, _ = curve_fit(modelGP, RedFreq, np.log10(Storage), p0=InitialGuess, absolute_sigma=True, sigma=Sigma)
        Coeff = [popt[0], popt[1], popt[2], popt[3]]
    # ------------------------------------------------------------------------------------------------------------------
    elif FitType == 2:      # Fit to the Storage modulus data (G').
        # Define the storage modulus model. 
        def modelGP(w, Gg, wc, m, k):
            G = modelGstar(w, Gg, wc, m, k)
            P = modelPhase(w, wc, m, k)
            return np.log10((10 ** G) * np.sin(np.deg2rad(P))) 
        # Use default initial Guess, if not provided. 
        if InitialGuess == []: 
            InitialGuess = [8.5, 2.0, 1.0, 0.4]
        # Using the "Curve fit" function, which uses unbound "Levenberg-Marquardt" algorithm and non-linear least 
        #   square fitting loss function.
        popt, _ = curve_fit(modelGP, RedFreq, np.log10(Loss), p0=InitialGuess, absolute_sigma=True, sigma=Sigma)
        Coeff = [popt[0], popt[1], popt[2], popt[3]]
    # ------------------------------------------------------------------------------------------------------------------
    elif FitType == 3:      # Fit to the G* and δ at the same time (7 parameters to fit).
        # Use default initial Guess, if not provided. 
        if InitialGuess == []: 
            InitialGuess = [8.5, 2.0, 1.0, 0.4, 2.0, 1.0, 0.4]
        # Using the "Curve fit" function, which uses unbound "Levenberg-Marquardt" algorithm and non-linear least 
        #   square fitting loss function.
        Gpopt, _ = curve_fit(modelGstar, RedFreq, np.log10(Gstar), p0=InitialGuess[:4], absolute_sigma=True,sigma=Sigma)
        Ppopt, _ = curve_fit(modelPhase, RedFreq, Phase, p0=InitialGuess[4:], absolute_sigma=True, sigma=Sigma)
        Coeff = [Gpopt[0], Gpopt[1], Gpopt[2], Gpopt[3], Ppopt[0], Ppopt[1], Ppopt[2]]
    # ------------------------------------------------------------------------------------------------------------------
    else: 
        raise Exception(f'Fit type is not recognized!! Please check the code manually!')
    # ------------------------------------------------------------------------------------------------------------------
    # Prepare outputs. 
    # Calculate the MAE loss function value. 
    Gstar_Pred = 10 ** modelGstar(RedFreq, Coeff[0], Coeff[1], Coeff[2], Coeff[3])
    Plot_Gstar = 10 ** modelGstar(Freq,    Coeff[0], Coeff[1], Coeff[2], Coeff[3])
    if len(Coeff) == 7:
        Phase_Pred = modelPhase(RedFreq, Coeff[4], Coeff[5], Coeff[6])
        Plot_Phase = modelPhase(Freq,    Coeff[4], Coeff[5], Coeff[6])
    else:
        Phase_Pred = modelPhase(RedFreq, Coeff[1], Coeff[2], Coeff[3])
        Plot_Phase = modelPhase(Freq,    Coeff[1], Coeff[2], Coeff[3])
    Phase_Pred[Phase_Pred > 90] = 90
    Plot_Phase[Plot_Phase > 90] = 90
    GP_Pred  = Gstar_Pred * np.cos(np.deg2rad(Phase_Pred))
    Plot_GP  = Plot_Gstar * np.cos(np.deg2rad(Plot_Phase))
    GPP_Pred = Gstar_Pred * np.sin(np.deg2rad(Phase_Pred))
    Plot_GPP = Plot_Gstar * np.sin(np.deg2rad(Plot_Phase))
    MLAE_Gstar = np.abs(np.log10(Gstar) - np.log10(Gstar_Pred)).sum() / NumData
    MLAE_Phase = np.abs(Phase - Phase_Pred).sum() / NumData
    MLAE_GP    = np.abs(np.log10(Storage) - np.log10(GP_Pred)).sum() / NumData
    MLAE_GPP   = np.abs(np.log10(Loss) - np.log10(GPP_Pred)).sum() / NumData
    # Extract the results and the model predictions. 
    Gg, wc, m, k = 10 ** Coeff[0], 10 ** Coeff[1], Coeff[2], Coeff[3]
    if len(Coeff) == 7:
        wcd, md, kd = Coeff[4], Coeff[5], Coeff[6]
    else:
        wcd, md, kd = wc, m, k
    # Calculating the indices for the fitted master curve. 
    Gwc = 10 ** modelGstar(wc, Coeff[0], Coeff[1], Coeff[2], Coeff[3])      # Crossover modulus (psi).
    R   = np.log10(Gg / Gwc)                                                # Rheological index. 
    # Store the results. 
    Res = {'MC_Coeff': [Gg, wc, m, k, wcd, md, kd],
           'MC_Coeff_Labels': ['Glassy Modulus (psi)', 'Crossover Frequency (rad/s)', 'Shape parameter, m', 
                               'Shape parameter, k', 'Crossover Frequency of δ MC (rad/s)', 'Shape parameter, md, δ MC', 
                               'Shape parameter, kd, δ MC'],
           'Plot_Freq': Freq, 'Plot_Gstar': Plot_Gstar, 'Plot_Phase': Plot_Phase, 
           'Plot_Storage': Plot_GP, 'Plot_Loss': Plot_GPP,
           'MC_Index': [R, wc, Gwc], 
           'MC_Index_Labels': ['Rheology Index', 'Crossover Frequency (rad/s)', 'Crossover Modulus (psi)'], 
           'MC_Error': [MLAE_Gstar, MLAE_Phase, MLAE_GP, MLAE_GPP], 
           'MC_Error_Labels': ['MLAE (|G*|)', 'MLAE (δ)', "MLAE (G')", 'MLAE (G")']}
    # Return the results. 
    return Res
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

