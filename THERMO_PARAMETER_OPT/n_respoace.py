import pickle
import os
import sys

# --- CONFIGURATION: PASTE YOUR FILE PATH HERE ---
PKL_FILE_PATH = "/data2/RANA/FINAL_RUN_7nc/FULL_RUN_1/RESULT_ANALYSIS_FOLDERS/response_surafce_analysis/ResponseSurfaces.pkl" 
# Index to inspect (e.g., 0 for the first ResponseSurface model)
MODEL_INDEX_TO_INSPECT = 0 
# --------------------------------------------------


def get_data_summary(value):
    """Generates a summary string for a value (e.g., type, shape, length)."""
    data_type = type(value).__name__
    summary = f"Type: {data_type}"
    
    # Check for common array/data structure types
    if hasattr(value, 'shape') and data_type == 'ndarray':
        summary += f", Shape: {value.shape}"
    elif isinstance(value, (list, tuple, set)):
        summary += f", Length: {len(value)}"

    # If it's a simple type or small, show a sample
    if not any(attr in summary for attr in ['Shape', 'Length']):
         # Show up to 80 characters of the string representation
         sample_value = str(value).replace('\n', ' ')
         summary += f", Value Sample: {sample_value[:80]}..."
         
    return summary


def inspect_single_object_attributes(obj, obj_name="Object"):
    """
    Inspects the internal attributes of a custom class instance or dictionary.
    """
    
    print(f"\n## Internal Inspection of '{obj_name}'")
    print("-" * 40)
    
    # Determine how to get the keys/attributes
    if isinstance(obj, dict):
        keys_to_inspect = obj.keys()
        data_source = obj
        structure_type = "Dictionary Key"
    elif hasattr(obj, '__dict__'):
        keys_to_inspect = obj.__dict__.keys()
        data_source = obj.__dict__
        structure_type = "Class Attribute"
    else:
        print(f"Cannot inspect internal structure: {get_data_summary(obj)}")
        return

    print(f"Found {len(keys_to_inspect)} {structure_type}s.")
    
    # Iterate through all keys/attributes
    for key in keys_to_inspect:
        try:
            # Retrieve the value
            if isinstance(data_source, dict):
                value = data_source[key]
            else:
                value = getattr(obj, key)
            
            key_repr = repr(key) 
            summary = get_data_summary(value)

            print(f"**{key_repr}**: {summary}")

        except Exception as e:
            print(f"  [ERROR] Could not inspect key {repr(key)}: {e}")


def inspect_pkl_content(file_path, model_index):
    """
    Loads a .pkl file, identifies the main structure, and then drills down 
    into a specific ResponseSurface model.
    """
    
    print("!!! SECURITY WARNING !!!")
    print("The 'pickle' module is not secure. Only load files from trusted sources.")
    print("-" * 50)
    
    if not os.path.exists(file_path):
        print(f"❌ Error: File not found at '{file_path}'. Please check the path.")
        return

    # 1. LOAD THE ENTIRE OBJECT
    try:
        with open(file_path, 'rb') as file:
            loaded_object = pickle.load(file)
            
        print(f"✅ Successfully loaded object from: {file_path}")
    except Exception as e:
        print(f"❌ Error loading file: {e}")
        return

    # 2. IDENTIFY MAIN STRUCTURE
    object_type = type(loaded_object)
    print("\n## Main Object Type and Structure")
    print("-" * 40)
    print(f"Main Object Type: {object_type.__name__}")
    
    if not isinstance(loaded_object, dict):
        print("Expected a dictionary structure based on previous output, but got a different type.")
        return
        
    print(f"This dictionary contains {len(loaded_object)} items (keys: 0 to {len(loaded_object) - 1}).")

    # 3. DRILL DOWN INTO A SINGLE RESPONSESURFACE OBJECT
    if model_index in loaded_object:
        target_model = loaded_object[model_index]
        
        # This calls the new function to inspect the internal attributes
        inspect_single_object_attributes(target_model, obj_name=f"ResponseSurface at Index {model_index}")
        
    else:
        print(f"\n❌ Error: Key '{model_index}' not found in the dictionary.")


if __name__ == "__main__":
    inspect_pkl_content(PKL_FILE_PATH, MODEL_INDEX_TO_INSPECT)
