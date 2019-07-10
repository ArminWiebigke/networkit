from plot_scripts.read_data import read_data
from plot_scripts.config import set_sns_style
from plot_scripts.create_plots import run

set_sns_style()

print("Reading results...")
data = read_data()

print("Creating plots...")
run(data)