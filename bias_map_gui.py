'''Code to produce bias maps using a GUI for inputs'''

#imports functions required
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import pandas as pd
from tqdm import tqdm
from tkinter import *
from tkinter import ttk


#defining various function required for the main run
def array_ln(v):
     
    '''Takes natural log of every item in an array'''
    
    ln_v = np.log(v)
    
    return ln_v
    
v_array_ln = np.vectorize(array_ln)



def model_pen(x, dma, dmb, sda, sdb):
    
    '''Calculates the model penetration curves'''
    
    a = 1 - norm.cdf(np.log(x), loc=np.log(dma), scale=np.log(sda))
    b = 1 - norm.cdf(np.log(x), loc=np.log(dmb), scale=np.log(sdb))
    penetration = a * b
    
    return penetration

v_model_pen = np.vectorize(model_pen)



def prob_den_calc_xl(dp_list, gsd, mmad):
    
    '''Calculates the probablility density for a given GSD and MMAD using the excel formula'''

    top = np.exp(-0.5 * np.square((np.log(dp_list / mmad)) / (np.log(gsd))))
    bottom = np.sqrt(6.28318) * np.log(gsd)
    # prob_den = top / (bottom * dp_list)
    prob_den = top / bottom
    
    return prob_den
    
v_prob_den_calc_xl = np.vectorize(prob_den_calc_xl)



def prob_den_calc(dp_list, gsd, mmad):
    
    '''Calculates the probablility density for a given GSD and MMAD using norm.pdf'''
    
    prob_den = norm.pdf(v_array_ln(dp_list), np.log(mmad), scale=np.log(gsd))
    
    return prob_den
    
v_prob_den_calc = np.vectorize(prob_den_calc)



def calc_resp(dp): # respirable as a fraction of inhalable
    
    '''Calculates the respirable fraction at differnt particle sizes'''

    frac = 1 - norm.cdf(np.log(dp), np.log(4.25), scale=np.log(1.5))

    return frac

v_calc_resp = np.vectorize(calc_resp)



def calc_inhal(dp): # inhalable
    
    '''Calculates the inhalible fraction at differnt particle sizes'''
    
    frac = 50 * (1 + np.exp(-0.06 * dp))
    
    return frac / 100

v_calc_inhal = np.vectorize(calc_inhal)



def bias_cal(prob_den, foam_pen, resp_pen):
    
    '''Calculates the bias for given probability and penatration distributions'''
    
    FMS_resp, FMS_foam = np.sum(prob_den * resp_pen), np.sum(prob_den * foam_pen)
    bias = ((FMS_foam - FMS_resp) / FMS_resp) * 100
    
    return bias



def relevant(compare, x, y, bias):
    
    '''Calculates if a point on the bias map is of interest'''
    
    relevant = False
    if (y >= (4 / 15) * x + (13 / 30)) and (y <= (-32 / 3) * x + (409 / 6) or y <= -14 * x + 76 or y <= -44 * x + 143.5):
        relevant = True
    if (compare == 'total' or compare == 'inhalible') and y >= 12 * x -8.5:
        relevant = False
        
    if relevant == True:
        return bias
    else:
        return None
    


def av_abs_bias(bias, bias_count, none_count):
    
    '''Prints out the average absolute bias in the area of intrest'''
    
    total_bias = 0
    for b_row in bias:
        for b in b_row:
            if b != None:
                total_bias += np.abs(b)
    av_bias = total_bias / (bias_count - none_count)
    to_print = "\nThe average absolute bias: " + str("{:.2f}".format(av_bias))
    L15 = Label(root, text=to_print).grid(row=11, column=0)
    
    return 



def percent_within_10(bias, bias_count, none_count):
    
    '''Prints out the percentage of the area of intrest within ten percent bias'''
    
    within_10 = 0
    for b_row in bias:
        for b in b_row:
            if b != None and np.abs(b) <= 10:
                within_10 += 1
    percent = (within_10 / (bias_count - none_count)) * 100
    to_print = "The percentage within 10% bias: " + str("{:.2f}".format(percent)) + "%"
    L16 = Label(root, text=to_print).grid(row=12, column=0, pady=10)
    
    return



def plot_pen(dp_list, foam_pen, inhalible, respirable, dp_min, dp_max):
    
    '''Plots the various penetration curves'''
    
    plt.figure(dpi=600)
    plt.plot(dp_list, foam_pen)
    plt.plot(dp_list, inhalible)
    plt.plot(dp_list, respirable)
    plt.legend(['foam', 'inhalible', 'respirable'])
    
    #Setting axes
    ax = plt.gca() 
    ax.set_xlim([dp_min, dp_max])
    ax.set_ylim([0.0, 1.0])
    plt.xlabel('Mass Median Aerodynamic Diameter / $\mu$m')
    plt.ylabel('Penetration')



def plot_bias(X, Y, Z, mmad_max, gsd_min, gsd_max):
    
    '''Created the contour plot with overlaid heatmap for the bias'''

    plt.figure(dpi=600)
    # print(X[20][7], Y[20][7], Z[20][7])
    CS1 = plt.contour(X, Y, Z) #plotting the contour map
    
    #Awkward bit of code to get round plotting area not of intrest as 'None' for contours but 0 for heatmap
    new_Z = []
    for row_z in Z:
        new_z = []
        for z in row_z:
            if z == None:
                z = 0
            new_z.append(z)
        new_Z.append(new_z)       
    plt.clabel(CS1)
    plt.imshow(new_Z, cmap='Blues_r', origin ='lower', aspect='auto', alpha=0.8, extent=(0, mmad_max, gsd_min, gsd_max))
    plt.colorbar()
    
    #Setting axes
    ax = plt.gca()
    ax.set_xlim([0, mmad_max])
    ax.set_ylim([gsd_min, gsd_max])
    plt.xlabel('Mass Median Aerodynamic Diameter / $\mu$m')
    plt.ylabel('Geometric Standard Deviation')
    
    #Saving and showing plot
    plt.savefig('Bias_Plot.png')
    plt.show()
    



def save_bias(params, Z):
        
    '''Saves the results of the bias caluclations and the parameters used to a csv file'''

    with open('Bias_Data.csv', 'w') as f:
        for key,value in params.items():
            f.write(str(key) + ',' + str(value) + '\n')
        
        f.write('\n')
        for i in Z:
            for j in i:
                if j == None:
                    f.write('###,')
                else:
                    f.write(str(j) + ",")
            f.write('\n')



def yes_no(yn):
    
    """Takes input of 'Yes' or 'No' and returns 'True' or 'False'"""
    
    if yn == "Yes":
        return True
    if yn == "No":
        return False



def myClick():
    
    """Gets the parameters from GUI and runs the calculation"""
    
    params = {
        'gsd_min': float(V1.get()),
        'gsd_max': float(V2.get()),
        'gsd_step': float(V3.get()),
        'mmad_min': float(V5.get()),
        'mmad_max': float(V6.get()),
        'mmad_step': float(V7.get()),
        'dp_min': float(V9.get()),
        'dp_max': float(V10.get()),
        'dp_step': float(V11.get()),
        'data_set': V12.get(), #Set as original or optimised
        'compare': V13.get(), #Set as total or repirable
        'set_dp': True, #Set as True to generate diameters or False to extract from Excel
        'area_of_intrest': yes_no(V14.get()), #Set as True to only see region of intrest or False to see all regions
        'max_value': 400 #Can be used to cap the bias map for ease of visualisation
    }
    
    main(**params)



def main(**args):
    
    '''Main code to cerate data needed to plot the bias map'''
    
    #Parameters being passed
    gsd_min, gsd_max, gsd_step = args['gsd_min'], args['gsd_max'], args['gsd_step']
    mmad_min, mmad_max, mmad_step = args['mmad_min'], args['mmad_max'], args['mmad_step']
    dp_min, dp_max, dp_step = args['dp_min'], args['dp_max'], args['dp_step']
    data_set, max_value = args['data_set'], args['max_value']
    compare, set_dp, area_of_intrest = args['compare'], args['set_dp'], args['area_of_intrest']
    
    #Creating the lists of values of mean diameter and geometric standard seviation to use
    gsd_range = np.arange(gsd_min, gsd_max + gsd_step, gsd_step)
    mmad_range = np.arange(mmad_min, mmad_max + mmad_step, mmad_step)
    
    #Creating the list of diameters to use
    if set_dp == True:
        dp_list = np.arange(dp_min, dp_max + dp_step, dp_step)
    #Or extracting the dp_list to use from the Excel spreadsheet
    else:
        workbook_name = 'HSE Metal foam penetration curves.xlsx'
        sheet_name = 'Tables'
        data = pd.read_excel(workbook_name, sheet_name=sheet_name).to_numpy()[5:,1:]
        dp_list = data[:,0]

    #Setting the parameters for the model foam penetration curves for either the original or optimised model curves
    if data_set == 'original':
        dma, dmb, sda, sdb = 7.30, 5.57, 2.02, 1.28
    elif data_set == 'optimised':
        dma, dmb, sda, sdb = 4.57, 5.37, 1.24, 1.76
    else:
        raise ValueError("data_set should be original or optimised")
    
    #Calculating the foam, inhalible and respirable penetration curves
    foam_pen = v_model_pen(dp_list, dma, dmb, sda, sdb)
    inhalible = v_calc_inhal(dp_list)
    respirable_of_inhalible = v_calc_resp(dp_list)
    respirable_of_total = inhalible * respirable_of_inhalible
    
    if compare == 'inhalible':
        compare_pen = respirable_of_inhalible
    elif compare == 'total':
        compare_pen = respirable_of_total
    else:
        raise ValueError("compare_pen should be inhalible or total")
    
    #Generating grid of values to use as inputs at each point on bias map
    X, Y = np.meshgrid(mmad_range, gsd_range)

    #Calculating the bias at each point on the map
    Z = [] 
    none_count = 0
    bias_count = 0
    #for gsd in tqdm(gsd_range):
    for gsd in gsd_range:
        z = []
        for mmad in mmad_range:
            prob_den = v_prob_den_calc(dp_list, gsd, mmad)
            bias = bias_cal(prob_den, foam_pen, compare_pen)
            bias_count +=1
            #Setting bias to 'None' if in an area not of interest
            if area_of_intrest == True:
                bias = relevant(compare, gsd, mmad, bias)
            #bias = np.log(bias)
            if bias == None:
                none_count +=1
            elif bias >= max_value:
                bias = max_value
            z.append(bias)
        Z.append(z)
        
    # Older bit of code for comparing the probability densities
    # prob_den = v_prob_den_calc(dp_list, 1.75, 7)
    # prob_den_xl = v_prob_den_calc_xl(dp_list, 1.75, 7)
    # plt.figure(dpi=600)
    # plt.plot(dp_list, prob_den)
    # plt.plot(dp_list, prob_den_xl)
    # plt.legend(['new method', 'excel method'])
    # plt.show()
    
    #Will plot the penetration curves or bias map as required
    #plot_pen(dp_list, foam_pen, inhalible, respirable, dp_min, dp_max)
    plot_bias(X, Y, Z, mmad_max, gsd_min, gsd_max)
    
    #Printing caulations around bias
    av_abs_bias(Z, bias_count, none_count)
    percent_within_10(Z, bias_count, none_count)
    
    #Saving the bias to csv file
    save_bias(args, Z)



if __name__ == "__main__":
    
    '''Setting up GUI'''
    
    root = Tk()
    root.title(" Bias Calculator")

    Heading = Label(root, text="Parameters for Calculation of Bias Map", font= ("TkDefaultFont", 10)).grid(row=0, column=0, pady=20)

    L0 = Label(root, text="Geometric Standard Deviation => ").grid(row=1, column=0)
    L1 = Label(root, text="Minimum:").grid(row=1, column=1)
    V1 = Entry(root, width=10)
    V1.grid(row=1, column=2)
    V1.insert(0, "1.75")
    L2 = Label(root, text="Maximum:").grid(row=1, column=3)
    V2 = Entry(root, width=10)
    V2.grid(row=1, column=4)
    V2.insert(0, "4.00")
    L3 = Label(root, text="Step:").grid(row=1, column=5)
    V3 = Entry(root, width=10)
    V3.grid(row=1, column=6)
    V3.insert(0, "0.05")

    L4 = Label(root, text="Mass Median Aerodynamic Diameter => ").grid(row=2, column=0)
    L5 = Label(root, text="Minimum:").grid(row=2, column=1)
    V5 = Entry(root, width=10)
    V5.grid(row=2, column=2)
    V5.insert(0, "0.5")
    L6 = Label(root, text="Maximum:").grid(row=2, column=3)
    V6 = Entry(root, width=10)
    V6.grid(row=2, column=4)
    V6.insert(0, "30.0")
    L7 = Label(root, text="Step:").grid(row=2, column=5)
    V7 = Entry(root, width=10)
    V7.grid(row=2, column=6)
    V7.insert(0, "0.5")

    L8 = Label(root, text="Calculation Particle Sizes (Resolution) => ").grid(row=3, column=0)
    L9 = Label(root, text="Minimum:").grid(row=3, column=1)
    V9 = Entry(root, width=10)
    V9.grid(row=3, column=2)
    V9.insert(0, "1.0")
    L10 = Label(root, text="Maximum:").grid(row=3, column=3)
    V10 = Entry(root, width=10)
    V10.grid(row=3, column=4)
    V10.insert(0, "18.0")
    L11 = Label(root, text="Step:").grid(row=3, column=5)
    V11 = Entry(root, width=10)
    V11.grid(row=3, column=6)
    V11.insert(0, "1.0")

    L12 = Label(root, text="Reference Curve (original/optimised)").grid(row=4, column=0)
    V12 = Entry(root, width=10)
    V12.grid(row=4, column=1, pady=10)
    V12.insert(0, "optimised")

    L13 = Label(root, text="Respirable as a fraction of (inhalible/total)").grid(row=5, column=0)
    V13 = Entry(root, width=10)
    V13.grid(row=5, column=1)
    V13.insert(0, "total")

    L14 = Label(root, text="Only Show Area of Intrest? (Yes/No)").grid(row=6, column=0)
    V14 = Entry(root, width=10)
    V14.grid(row=6, column=1, pady=10)
    V14.insert(0, "Yes")
    
    myButton = Button(root, text="Run Calulation", bg="#33b249", command=myClick)
    myButton.grid(row=10, column=0, pady=20)

    root.mainloop()





