import numpy as np

def get_rxn_list(sens_data, unsrt_data):
    rxn_dict = {}
    sc_dict = {}
    for case in sens_data:
        rxn_list =[]
        sc_list = []
        for rxn in sens_data[case]:
            rxn_list.append(rxn)
            sc_list.append(sens_data[case][rxn])
        rxn_dict[case] = rxn_list
        sc_dict[case] = sc_list
        #print(sc_list)
    new_rxn_dict = {}
    new_sc_dict = {}
    for case in rxn_dict:
        new_rxn_list = []
        new_sc_list = []
        rxn_list = rxn_dict[case]
        _rxn_list = [rxn for rxn in unsrt_data]
        for index,rxn in enumerate(rxn_list):
            for _rxn in _rxn_list:
                #print(rxn,_rxn)
                if _rxn.split(":")[0] == rxn:
                     new_rxn_list.append(_rxn)
                     new_sc_list.append(sc_dict[case][index])
        new_rxn_dict[case] = new_rxn_list
        new_sc_dict[case] = new_sc_list
        #print(new_sc_list)    
    return new_rxn_dict, new_sc_dict

