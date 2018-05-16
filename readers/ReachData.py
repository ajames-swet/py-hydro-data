import pandas as pd
import math
import numpy as np

"""Useful info: http://earthpy.org/pandas-basics.html"""


class WAMReachData:
    """
    WAM Reach data class: contains data to manipulate reach data from WAM simulations. 
    Upon reading the file, sums various forms, e.g., Soluble N is NO3 + Organic N + dissolved NH3
    """
    def __init__(self, reach_data = None, filename = None):
        self.reach_data = reach_data
        self.filename = filename
        if (self.filename != None):
            self.read_reach(self.filename)


    def read_reach(self, filename):
        """
        Reach a WAM reach from the supplied filename
        """
        self.filename = filename
        df = pd.read_csv(self.filename, parse_dates=True)
        # Create date-time index
        # See https://chrisalbon.com/python/pandas_time_series_basics.html
        df['SimDate'] = pd.to_datetime(df['SimDate'])
        df.index = df['SimDate']
        del df['SimDate']
        daily_vol = 24*3600*df.Flow
        df['TotalSolN'] = df.SolNO3Out + df.SolNH4Out + df.SolOrgNOut
        df['TotalPartN'] = df.SedNH4Out + df.SedOrgNOut
        df['TotalNH4'] = df.SolNH4Out + df.SedNH4Out
        df['TotalN'] = df.TotalSolN + df.TotalPartN
        df['TotalKN'] = df.SolNH4Out + df.SolOrgNOut + df.SedNH4Out + df.SedOrgNOut
        df['TotalP'] = df.SolPOut + df.SedPOut
        # Loads are in kg/d
        df['SolNO3Load'] = daily_vol*df.SolNO3Out/1000.0
        df['SolNH4Load'] = daily_vol*df.SolNH4Out/1000.0
        df['SolOrgNLoad'] = daily_vol*df.SolOrgNOut/1000.0
        df['SolPLoad'] = daily_vol*df.SolPOut/1000.0
        df['SedNH4Load'] = daily_vol*df.SedNH4Out/1000.0
        df['SedOrgNLoad'] = daily_vol*df.SedOrgNOut/1000.0
        df['SedPLoad'] = daily_vol*df.SedPOut/1000.0
        df['TotalSolNLoad'] = daily_vol*df.TotalSolN/1000.0
        df['TotalPartNLoad'] = daily_vol*df.TotalPartN/1000.0
        df['TotalNLoad'] = daily_vol*df.TotalN/1000.0
        df['TotalPLoad'] = daily_vol*df.TotalP/1000.0
        df['BODLoad'] = daily_vol*df.BODOut/1000.0
        df['TSSLoad'] = daily_vol*df.TSSOut/1000.0

        df = df.drop(df.index[[0]])
        self.reach_data = df

    def calculate_flow_weighted_conc(self):
        total_vol = 24.0*3600*(self.reach_data.Flow.sum())
        sol_no3 = 1000.0*(self.reach_data.SolNO3Load.sum())/total_vol
        sol_nh4 = 1000.0*(self.reach_data.SolNH4Load.sum())/total_vol
        sol_org_n = 1000.0*(self.reach_data.SolOrgNLoad.sum())/total_vol
        sol_p = 1000.0*(self.reach_data.SolPLoad.sum())/total_vol
        sed_nh4 = 1000.0*(self.reach_data.SedNH4Load.sum())/total_vol
        sed_org_n = 1000.0*(self.reach_data.SedOrgNLoad.sum())/total_vol
        sed_p = 1000.0*(self.reach_data.SedPLoad.sum())/total_vol
        tot_n = 1000.0*(self.reach_data.TotalNLoad.sum())/total_vol
        tot_p = 1000.0*(self.reach_data.TotalPLoad.sum())/total_vol
        d = {'Total Vol': [total_vol], 'SolNO3': [sol_no3], 'SolNH4': [sol_nh4], 'SolOrgN': [sol_org_n], 'SolP': [sol_p], 
        'SedNH4': [sed_nh4], 'SedOrgN': [sed_org_n], 'SedP': [sed_p], 'TotalN': [tot_n], 'TotalP': [tot_p]}
        df = pd.DataFrame(data = d)
        return df


    def get_daily_values(self, col):
        return (self.reach_data[col]).resample('D').mean()

    def calculate_loads(self):
        daily_vol = 24*3600*self.reach_data.Flow
        self.reach_data['TotalSolN'] = self.reach_data.SolNO3Out + self.reach_data.SolNH4Out + self.reach_data.SolOrgNOut
        self.reach_data['TotalPartN'] = self.reach_data.SedNH4Out + self.reach_data.SedOrgNOut
        self.reach_data['TotalN'] = self.reach_data.TotalSolN + self.reach_data.TotalPartN
        self.reach_data['TotalP'] = self.reach_data.SolPOut + self.reach_data.SedPOut
        # Loads are in kg/d
        self.reach_data['SolNO3Load'] = daily_vol*self.reach_data.SolNO3Out/1000.0
        self.reach_data['SolNH4Load'] = daily_vol*self.reach_data.SolNH4Out/1000.0
        self.reach_data['SolOrgNLoad'] = daily_vol*self.reach_data.SolOrgNOut/1000.0
        self.reach_data['SolPLoad'] = daily_vol*self.reach_data.SolPOut/1000.0
        self.reach_data['SedNH4Load'] = daily_vol*self.reach_data.SedNH4Out/1000.0
        self.reach_data['SedOrgNLoad'] = daily_vol*self.reach_data.SedOrgNOut/1000.0
        self.reach_data['SedPLoad'] = daily_vol*self.reach_data.SedPOut/1000.0
        self.reach_data['TotalSolNLoad'] = daily_vol*self.reach_data.TotalSolN/1000.0
        self.reach_data['TotalPartNLoad'] = daily_vol*self.reach_data.TotalPartN/1000.0
        self.reach_data['TotalNLoad'] = daily_vol*self.reach_data.TotalN/1000.0
        self.reach_data['TotalPLoad'] = daily_vol*self.reach_data.TotalP/1000.0
        self.reach_data['BODLoad'] = daily_vol*self.reach_data.BODOut/1000.0
        self.reach_data['TSSLoad'] = daily_vol*self.reach_data.TSSOut/1000.0


    def calculate_total_concentrations(self):
        self.reach_data['TotalSolN'] = self.reach_data.SolNO3Out + self.reach_data.SolNH4Out + self.reach_data.SolOrgNOut
        self.reach_data['TotalPartN'] = self.reach_data.SedNH4Out + self.reach_data.SedOrgNOut
        self.reach_data['TotalN'] = self.reach_data.TotalSolN + self.reach_data.TotalPartN
        self.reach_data['TotalP'] = self.reach_data.SolPOut + self.reach_data.SedPOut

    def recalculate_concentrations(self):
#        flow_data = self.reach_data.Flow[not np.isclose(self.reach_data.Flow, 0.0).all()]
        flow_data = self.reach_data.Flow[abs(self.reach_data.Flow) > 0.0]

        daily_vol = 24*3600*flow_data
        # print("Daily vol:")
        # print(daily_vol.tail())
        # Multiply daily loads by 1000 (kg -> g) and divide by vol in m^3
        self.reach_data.SolNO3Out = self.reach_data['SolNO3Load']*1000.0/daily_vol
        self.reach_data.SolNH4Out = self.reach_data['SolNH4Load']*1000.0/daily_vol
        self.reach_data.SolOrgNOut = self.reach_data['SolOrgNLoad']*1000.0/daily_vol
        self.reach_data.SolPOut = self.reach_data['SolPLoad']*1000.0/daily_vol
        self.reach_data.SedNH4Out = self.reach_data['SedNH4Load']*1000.0/daily_vol
        self.reach_data.SedOrgNOut = self.reach_data['SedOrgNLoad']*1000.0/daily_vol
        self.reach_data.SedPOut = self.reach_data['SedPLoad']*1000.0/daily_vol
        self.reach_data.BODOut = self.reach_data['BODLoad']*1000.0/daily_vol
        self.reach_data.TSSOut = self.reach_data['TSSLoad']*1000.0/daily_vol
        # print("SOl org N final!")
        # print(self.reach_data.SolOrgNOut.head())
        self.calculate_total_concentrations()

    # See https://chrisalbon.com/python/pandas_group_data_by_time.html
    def get_monthly_sums(self, col, conversion = 1.0):
        return conversion*(self.reach_data[col]).resample('M').sum()

    def get_max_monthly_value(self, col, conversion = 1.0):
        return conversion*(self.reach_data[col]).resample('M').max()

    def get_annual_sums(self, col, conversion=1.0):
        return conversion * (self.reach_data[col]).resample('A').sum()
        # Q: does 'D' resample by day? So I can accumulate multiple grabs on a single day, and then pd.merge()?

    def get_difference(self, col, other_dataframe):
        return self.reach_data[col] - other_dataframe.reach_data[col]

    def merge_reach(self, next_reach):
        # Does this copy or reference?
        merged_data = WAMReachData(reach_data=self.reach_data.copy(), filename=None)
        merged_data.reach_data.Flow += next_reach.reach_data.Flow
        # Average concentrations, but recalculate flow-weighted later (if total flow is not zero)
        merged_data.reach_data.SolNO3Out = 0.5*(merged_data.reach_data.SolNO3Out + next_reach.reach_data.SolNO3Out)
        merged_data.reach_data.SolNH4Out = 0.5*(merged_data.reach_data.SolNH4Out + next_reach.reach_data.SolNH4Out)
        merged_data.reach_data.SolOrgNOut = 0.5*(merged_data.reach_data.SolOrgNOut + next_reach.reach_data.SolOrgNOut)
        merged_data.reach_data.SolPOut = 0.5*(merged_data.reach_data.SolPOut + next_reach.reach_data.SolPOut)
        merged_data.reach_data.SedNH4Out = 0.5*(merged_data.reach_data.SedNH4Out + next_reach.reach_data.SedNH4Out)
        merged_data.reach_data.SedOrgNOut = 0.5*(merged_data.reach_data.SedOrgNOut + next_reach.reach_data.SedOrgNOut)
        merged_data.reach_data.SedPOut = 0.5*(merged_data.reach_data.SedPOut + next_reach.reach_data.SedPOut)
        merged_data.reach_data.BODOut = 0.5*(merged_data.reach_data.BODOut + next_reach.reach_data.BODOut)
        merged_data.reach_data.TSSOut = 0.5*(merged_data.reach_data.TSSOut + next_reach.reach_data.TSSOut)
        merged_data.reach_data['SolNO3Load'] += next_reach.reach_data['SolNO3Load']
        merged_data.reach_data['SolNH4Load'] += next_reach.reach_data['SolNH4Load']
        merged_data.reach_data['SolOrgNLoad'] += next_reach.reach_data['SolOrgNLoad']
        merged_data.reach_data['SolPLoad'] += next_reach.reach_data['SolPLoad']
        merged_data.reach_data['SedNH4Load'] += next_reach.reach_data['SedNH4Load']
        merged_data.reach_data['SedOrgNLoad'] += next_reach.reach_data['SedOrgNLoad']
        merged_data.reach_data['SedPLoad'] += next_reach.reach_data['SedPLoad']
        merged_data.reach_data['TotalSolNLoad'] += next_reach.reach_data['TotalSolNLoad']
        merged_data.reach_data['TotalPartNLoad'] += next_reach.reach_data['TotalPartNLoad']
        merged_data.reach_data['TotalNLoad'] += next_reach.reach_data['TotalNLoad']
        merged_data.reach_data['TotalPLoad'] += next_reach.reach_data['TotalPLoad']
        merged_data.reach_data['BODLoad'] += next_reach.reach_data['BODLoad']
        merged_data.reach_data['TSSLoad'] += next_reach.reach_data['TSSLoad']
        merged_data.recalculate_concentrations()
        return merged_data

    def merge_reaches(self, reach_list):
        # Does this copy or reference?
#        merged_data = ReachData(reach_data=self.reach_data.copy(), filename=None)
        merged_data = WAMReachData(reach_data=None, filename=self.filename)
        num_reaches = 1.0 + float(len(reach_list))
        for next_reach in reach_list:
            merged_data.reach_data.Flow += next_reach.reach_data.Flow
            merged_data.reach_data.SolNO3Out += next_reach.reach_data.SolNO3Out
            merged_data.reach_data.SolNH4Out += next_reach.reach_data.SolNH4Out
            merged_data.reach_data.SolOrgNOut += next_reach.reach_data.SolOrgNOut
            merged_data.reach_data.SolPOut += next_reach.reach_data.SolPOut
            merged_data.reach_data.SedNH4Out += next_reach.reach_data.SedNH4Out
            merged_data.reach_data.SedOrgNOut += next_reach.reach_data.SedOrgNOut
            merged_data.reach_data.SedPOut += next_reach.reach_data.SedPOut
            merged_data.reach_data.BODOut += next_reach.reach_data.BODOut
            merged_data.reach_data.TSSOut += next_reach.reach_data.TSSOut
            merged_data.reach_data['SolNO3Load'] += next_reach.reach_data['SolNO3Load']
            merged_data.reach_data['SolNH4Load'] += next_reach.reach_data['SolNH4Load']
            merged_data.reach_data['SolOrgNLoad'] += next_reach.reach_data['SolOrgNLoad']
            merged_data.reach_data['SolPLoad'] += next_reach.reach_data['SolPLoad']
            merged_data.reach_data['SedNH4Load'] += next_reach.reach_data['SedNH4Load']
            merged_data.reach_data['SedOrgNLoad'] += next_reach.reach_data['SedOrgNLoad']
            merged_data.reach_data['SedPLoad'] += next_reach.reach_data['SedPLoad']
            merged_data.reach_data['TotalSolNLoad'] += next_reach.reach_data['TotalSolNLoad']
            merged_data.reach_data['TotalPartNLoad'] += next_reach.reach_data['TotalPartNLoad']
            merged_data.reach_data['TotalNLoad'] += next_reach.reach_data['TotalNLoad']
            merged_data.reach_data['TotalPLoad'] += next_reach.reach_data['TotalPLoad']
            merged_data.reach_data['BODLoad'] += next_reach.reach_data['BODLoad']
            merged_data.reach_data['TSSLoad'] += next_reach.reach_data['TSSLoad']

        # Average concentrations, but recalculate flow-weighted later (if total flow is not zero)
        #print ("Temp Flow final:")
        #print(merged_data.reach_data.Flow.tail())
        #print ("Temp org n load final:")
        #print(merged_data.reach_data['SolOrgNLoad'].tail())
        merged_data.reach_data.SolNO3Out /= num_reaches
        merged_data.reach_data.SolNH4Out /= num_reaches
        merged_data.reach_data.SolOrgNOut /= num_reaches
        merged_data.reach_data.SolPOut /= num_reaches
        merged_data.reach_data.SedNH4Out /= num_reaches
        merged_data.reach_data.SedOrgNOut /= num_reaches
        merged_data.reach_data.SedPOut /= num_reaches
        merged_data.reach_data.BODOut /= num_reaches
        merged_data.reach_data.TSSOut /= num_reaches
        merged_data.recalculate_concentrations()
        return merged_data

    def save_to_csv(self, filename):
        self.reach_data.to_csv(filename)

    def merge_nps_files(self, nps_data):
        merged_data = WAMReachData(reach_data=None, filename=self.filename)
        day_to_sec = 1.0/24.0/3600.0
        g_day_to_kg_sec = 1.0/24.0/3600.0/1000.0

        merged_data.reach_data.Flow += day_to_sec*(nps_data.reach_nps_data['Flow_Surf'] + nps_data.reach_nps_data['Flow_GW'])
        no3_load = g_day_to_kg_sec*(nps_data.reach_nps_data['SolNO3_Surf'] + nps_data.reach_nps_data['SolNO3_GW'])
        nh4_sol_load = g_day_to_kg_sec*(nps_data.reach_nps_data['SolNH4_Surf'] + nps_data.reach_nps_data['SolNH4_GW'])
        nh4_part_load = g_day_to_kg_sec*(nps_data.reach_nps_data['SedNH4_Surf'])
        org_n_sol_load = g_day_to_kg_sec*(nps_data.reach_nps_data['SolOrgN_Surf'] + nps_data.reach_nps_data['SolOrgN_GW'])
        org_n_part_load = g_day_to_kg_sec*(nps_data.reach_nps_data['SedOrgN_Surf'])
        p_sol_load = g_day_to_kg_sec*(nps_data.reach_nps_data['SolP_Surf'] + nps_data.reach_nps_data['SolP_GW'])
        p_sed_load = g_day_to_kg_sec*(nps_data.reach_nps_data['SedP_Surf'])

        merged_data.reach_data['SolNO3Load'] += no3_load
        merged_data.reach_data['SolNH4Load'] += nh4_sol_load
        merged_data.reach_data['SolOrgNLoad'] += org_n_sol_load
        merged_data.reach_data['SolPLoad'] += p_sol_load
        merged_data.reach_data['SedNH4Load'] += nh4_part_load
        merged_data.reach_data['SedOrgNLoad'] += org_n_part_load
        merged_data.reach_data['SedPLoad'] += p_sed_load
        merged_data.reach_data['BODLoad'] += g_day_to_kg_sec*(nps_data.reach_nps_data['BOD_Surf'])
        merged_data.reach_data['TSSLoad'] += g_day_to_kg_sec*(nps_data.reach_nps_data['TSS_Surf'])
        merged_data.reach_data['TotalSolNLoad'] += no3_load + nh4_sol_load + org_n_sol_load
        merged_data.reach_data['TotalPartNLoad'] += nh4_part_load + org_n_part_load
        merged_data.reach_data['TotalNLoad'] += no3_load + nh4_sol_load + org_n_sol_load + nh4_part_load + org_n_part_load
        merged_data.reach_data['TotalPLoad'] += p_sol_load + p_sed_load
        merged_data.recalculate_concentrations()
        return merged_data

class WAMReachNPSData:
    nps_names = ['Date', 'Flow_Surf', 'TSS_Surf', 'SolNO3_Surf', 'SolNH4_Surf', 'SolOrgN_Surf', 'SolP_Surf',
                 'SedNH4_Surf', 'SedOrgN_Surf', 'SedP_Surf', 'BOD_Surf', 'Flow_GW', 'SolNO3_GW', 'SolNH4_GW',
                 'SolOrgN_GW', 'SolP_GW']

    """Reach data class: contains data to manipulate reach data"""
    def __init__(self, reach_data = None, filename = None):
        self.reach_nps_data = reach_data
        self.filename = filename
        if (self.filename != None):
            self.read_reach(self.filename)

    # Default header:
    # Date,Flow_Surf(m3/d),TSS(g/d)_Surf,SolNO3_Surf(g/d),SolNH4_Surf(g/d),SolOrgN_Surf(g/d),SolP_Surf(g/d),SedNH4_Surf(g/d),SedOrgN_Surf(g/d),SedP_Surf(g/d),BOD_Surf(g/d),Flow_GW(m3/d),SolNO3_GW(g/d),SolNH4_GW(g/d),SolOrgN_GW(g/d),SolP_GW(g/d)
    def read_reach(self, filename):
        self.filename = filename
        parser = lambda date: pd.datetime.strptime(date, '%m/%d/%Y')
        df = pd.read_csv(self.filename, parse_dates=True, header=0, names=self.nps_names)
        # Create date-time index
        # See https://chrisalbon.com/python/pandas_time_series_basics.html
        df['Date'] = pd.to_datetime(df['Date'])
        df.index = df['Date']
        del df['Date']
        self.reach_nps_data = df

    def get_daily_values(self, col):
        return (self.reach_nps_data[col]).resample('D').mean()


class WQSingleComponent:
    """Reach data class: contains data to manipulate reach data"""
    wq_names = ['Station', 'Date', 'WQ_Value', 'DateBin', 'Unknown']

    def __init__(self, wq_data=None, filename=None):
        self.wq_data = wq_data
        self.filename = filename
        if (self.filename != None):
            self.read_data(self.filename)

    def read_data(self, filename):
        self.filename = filename
        df = pd.read_csv(self.filename, parse_dates=True, header=None, names=self.wq_names)

        # Create date-time index
        # See https://chrisalbon.com/python/pandas_time_series_basics.html
        df['Date'] = pd.to_datetime(df['Date'])
        df.index = df['Date']
        del df['Date']
        self.wq_data = df
