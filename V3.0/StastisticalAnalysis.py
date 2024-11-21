import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from fpdf import FPDF

class StatisticalAnalysis:
    def __init__(self, data,case_index):
        self.data = data
        self.case_index = case_index

    def z_score(self):
        return stats.zscore(self.data)

    def auto_split_groups(self):
        group1 = [self.data[i] for i in range(len(self.data)) if i % 2 == 0]
        group2 = [self.data[i] for i in range(len(self.data)) if i % 2 != 0]
        return group1, group2

    def t_test(self):
        group1, group2 = self.auto_split_groups()
        return stats.ttest_ind(group1, group2)

    def iqr(self):
        Q1, Q3 = np.percentile(self.data, [25, 75])
        return Q3 - Q1, Q1, Q3

    def extreme_values(self):
        IQR, Q1, Q3 = self.iqr()
        lower_bound, upper_bound = Q1 - 1.5 * IQR, Q3 + 1.5 * IQR
        return [x for x in self.data if x < lower_bound or x > upper_bound]

    def plot_data_spread(self, save_path):
        # Plot boxplot and save as image
        plt.figure(figsize=(8, 6))
        plt.boxplot(self.data, vert=False, patch_artist=True, boxprops=dict(facecolor='lightblue'))
        plt.title("Data Spread Visualization (Boxplot)")
        plt.xlabel("Values")
        plt.savefig(save_path)  # Save the plot to a file
        plt.close()

    def plot_parity(self, save_path):
        outliers, lower_bound, upper_bound = self.extreme_values()
        
        # Plotting
        plt.figure(figsize=(8, 6))
        plt.scatter(range(len(self.data)), self.data, label='Data Points', color='blue')
        
        # Highlight outliers
        outlier_indices = [i for i, val in enumerate(self.data) if val < lower_bound or val > upper_bound]
        plt.scatter(outlier_indices, [self.data[i] for i in outlier_indices], color='red', label='Outliers', s=100)

        plt.axhline(y=lower_bound, color='orange', linestyle='--', label='Lower Bound')
        plt.axhline(y=upper_bound, color='green', linestyle='--', label='Upper Bound')
        plt.title("Parity Plot with Outliers")
        plt.xlabel("Index")
        plt.ylabel("Values")
        plt.legend()
        plt.savefig(save_path)  # Save the plot to a file
        plt.close()
    
    def generate_pdf(self, save_path):
        # Create a PDF document
        pdf = FPDF()
        pdf.set_auto_page_break(auto=True, margin=15)
        pdf.add_page()

        # Add title
        pdf.set_font("Arial", "B", 16)
        pdf.cell(200, 10, "Statistical Analysis Report", ln=True, align="C")

        # Add Z-scores
        z_scores = self.z_score()
        pdf.set_font("Arial", "", 12)
        pdf.cell(200, 10, "Z-scores:", ln=True)
        pdf.multi_cell(0, 10, ', '.join([f"{z:.2f}" for z in z_scores]))

        # Add T-test
        t_stat, p_value = self.t_test()
        pdf.cell(200, 10, "T-test Results:", ln=True)
        pdf.multi_cell(0, 10, f"T-statistic: {t_stat:.4f}, P-value: {p_value:.4f}")

        # Add IQR
        iqr_value, Q1, Q3 = self.iqr()
        pdf.cell(200, 10, "Interquartile Range (IQR):", ln=True)
        pdf.multi_cell(0, 10, f"IQR: {iqr_value:.2f}, Q1: {Q1:.2f}, Q3: {Q3:.2f}")

        # Add Extreme Values
        outliers = self.extreme_values()
        pdf.cell(200, 10, "Extreme Values (Outliers):", ln=True)
        pdf.multi_cell(0, 10, ', '.join(map(str, outliers)))


        # Add plot
        pdf.cell(200, 10, "Data Spread (Boxplot):", ln=True)
        plot_path = f"../Plots/stats_spread_{self.case_index}.png"  # File path for saving plot
        self.plot_data_spread(plot_path)  # Create and save the plot
        pdf.image(plot_path, x=10, y=pdf.get_y(), w=180)
	
	#Parity plot
        pdf.cell(200, 10, "Parity Plot with Outliers:", ln=True)
        plot_path_parity = "../Plots/Parity_plot_with_outliers_{self.case_index}.png"
        self.plot_parity(plot_path_parity)  # Create and save the parity plot
        pdf.image(plot_path_parity, x=10, y=pdf.get_y(), w=180)
        # Save the PDF
        pdf.output(save_path)
        print(f"PDF report saved to {save_path}")

