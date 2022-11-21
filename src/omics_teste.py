import omicscope as omics

from copy import copy
#general = Omicscope(Table = 'C:/Users/Guilherme/Desktop/general.xls', ControlGroup = None, Method = 'General')
progenesis = omics.Omicscope('C:/Users/Guilherme/Desktop/progenesis.csv',ControlGroup = None, Method = 'Progenesis')
progenesis
def expression(self):
        """Joins the technical replicates and organizes biological
        conditions.
        """
        from copy import copy
        pdata = []
        for i in self.pdata.columns:
            pdata.append(self.pdata[i])
        expression = self.assay.T
        expression = expression.set_index(pdata).T

        rdata = []
        for i in self.rdata.columns:
            rdata.append(self.rdata[i])
        expression = expression.set_index(rdata)


        # Melt expression data to get all values for each sample and each gene
        expression = expression.melt(ignore_index=False).reset_index()
        # Getting the mean of technical replicates for each biological replicate
        var = list(self.rdata.columns)
        pdata_columns = list(self.pdata.columns)
        pdata_columns = list(set(pdata_columns) - set(['Sample', 'Samples', 'TechRep']))
        for i in pdata_columns:
            var.append(i)

        expression = expression.groupby(var).mean(numeric_only= True).reset_index().rename(
            columns={'value': 'abundance'})

        expression['Sample'] = 'BioRep' + \
            expression['Biological'].astype(str) + \
            '.' + expression['Condition']
        
        # New pdata
        pdata_columns.append('Sample')
        pdata = expression[pdata_columns]
        self.new_pdata = pdata

        # New rdata 
        rdata = expression[copy(self.rdata.columns)]
        self.new_rdata = rdata
        #  data with just biological replicates abundance
        expression = expression[['Sample', 'Accession', 'abundance']]
        expression = expression.pivot(columns='Sample', index='Accession', values='abundance')
        return(expression)

#result = expression(copy(general))
#result