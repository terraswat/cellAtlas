"""
Determine the data types of cell/sample annotation data.

Expects a pandas data frame where the columns are attributes and the
rows are cell/sample ids.

Some test data

attribute_df = pd.DataFrame(
    [["a", 1, 5.0, "yes"], ["b",2,3,"yes"], ["c",1,6, "no"]],
    columns=["no", "yep", "nope", "yes"]
)
should produce 
{'bin': ['yep', 'yes'], 'cat': ['no'], 'cont': ['nope']}

"""

def data_type_dictionary(attribute_df):
    data_type_dict = {}
    two_values = lambda x: len(x.dropna().unique()) == 2
    # True/false series for which attributes are binary.
    bin_mask = attribute_df.apply(two_values, axis=0).values
    
    # Grab the binary attributes.
    bin_attributes = attribute_df.columns[bin_mask].tolist()
    # Grab the rest of the attibutes.
    attribute_df = attribute_df.iloc[:,~bin_mask]

    # Categorical are strings and will be read in as type
    # object in pandas.
    cat_mask = attribute_df.dtypes == object
    
    cat_attributes = attribute_df.columns[cat_mask].tolist()

    # The rest are continuous attributes.
    cont_attributes = attribute_df.columns[~cat_mask].tolist()
    
    data_type_dict = {
        "bin" : bin_attributes,
        "cat" : cat_attributes,
        "cont" : cont_attributes
    }
    return data_type_dict 

def main():
    return 0
