#Import sys for obtaining command line args.
import numpy as np
import pickle as pkl

#Import matplotlib for plots.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sys.path.append(os.path.abspath("../common_scripts"))
import region_defs
import wig_and_signal_utils as wsu

"""
The following function retrieves all signals of a specified RE
in the given pickle file.
"""
def get_all_signals_for_re(shapes_file, re, promoter, enhancer, repressor, weak, is_peas):
    
    # Read the file.
    shapes = pkl.load(open(shapes_file, "rb"))
    signals = []
    positions = []
   
    # Build a list of signals for the RE.
    for shape in shapes:
        annotation = wsu.get_annotation(shape.signals, promoter, enhancer, repressor, weak, is_peas):
        if annotation == re:
            for i in range(len(shape.signals)):
                positions.append(i)
                signals.append(shape.signals[i])
    
    all_signals = pd.DataFrame({"Signal": signals, "Position": positions}, columns = ["Position", "Signal"])
    return all_signals
    

"""
Plot the average signal shape for the given data with
a confidence interval.
"""
def make_one_plot(title, signals):

    sns.lineplot(x="Position", y="Signal", data=signals)
    plt.title(title)

if __name__ == "__main__":
    main()