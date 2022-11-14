def generate_initiation_interface():
    # generates a GUI with an 
    # IntText inut for the STudent_ID,
    # Button to start iniate_Params_HW(#Homeworknumber)() with
    # output to print the values from iniate_Params_HW(#Homeworknumber)() in
    StudID = widgets.IntText( description='Student ID:', value=0 )
    button_StudID   = widgets.Button( description='Generate Parameters' )
    output = widgets.Output()
    
    global vb
    vb=widgets.VBox([StudID,button_StudID, output])
    return

##########
##########
##########
# functions for compute_Params() for every homework_notebook
##########
##########
##########

def compute_Params_HW3(b):
    # diable Button from initatl GUI
    # create hw3 in global namespace
    # store homework_three instance in hw3
    # print in output widget of inital GUI
    vb.children[1].disabled = True
    global hw3
    hw3 = homework_three(vb.children[0].value)
    with vb.children[2]:
        print("##########\n\nParameters for Question 1:\n\n"\
        "Your glucose yield Y_(XS) is {} C-mol glucose/C-mol biomass\n\n"\
        "Your biomass composition is CH_({})O_({})N_({})\n\n##########".format(hw3.newYxs,hw3.newH,hw3.newO,hw3.newN))
        print("##########\n\nParameters for Question 2:\n\n"\
        "Your O_2 Exchange-Rate is {} mol/h\n\n"\
        "Your CO_2 Exchange-Rate is {} mol/h\n\n##########".format(hw3.newRO2,hw3.newRCO2))

#####
#####
#####

def compute_Params_HW4(b):
    # diable Button from initatl GUI
    # create hw4 in global namespace
    # store homework_three instance in hw3
    # print in output widget of inital GUI
    vb.children[1].disabled = True
    global hw4
    hw4 = homework_four(vb.children[0].value)
    with vb.children[2]:
        #NOTE: Chemical Compositions are, as of yet, not actually modified in HW4 
        #      but are treated like it for the sake of the students
        print("##########\n\nParameters for Question 1:\n\n"\
        "Approximate Chemical Compostion of Fats: CH_(1.92)O_(0.12)\n" \
        "Approximate Chemical Compostion of Proteins: CH_(1.57)O_(0.32)N_(0.26)\n\n" \
        "Heat of combustion for protein: {} kcal/g,\n" \
        "Heat of combustion for Carbohydrates: {} kcal/g,\n" \
        "Heat of combustion for Fat: {} kcal/g\n\n##########".format(hw4.newCoP,hw4.newCoC,hw4.newCoF))
        print("##########\n\nParameters for Question 3:\n\n"\
        "Oxygen consumption rate: {} mol/h,\n" \
        "Carbon dioxide production rate: {} mol/h,\n" \
        "Nitrogen secretion rate: {} g N/h\n\n##########".format(hw4.newOCR,hw4.newCPR, hw4.newNSR))