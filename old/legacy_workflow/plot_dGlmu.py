import ROOT



def plot_comparison(Data_file, MC_file, plot1, plot2, tittle, x_axis_title, y_axis_title,y_low_limit,y_max_limit, file_name):
    

    # Open the two ROOT files
    file1_data = ROOT.TFile.Open(Data_file, "READ")
    file2_mc = ROOT.TFile.Open(MC_file, "READ")
    # Retrieve graphs by name Cosmics_muons_MC_DGL.root
    graph1_data = file1_data.Get(plot1)
    graph2_mc = file2_mc.Get(plot2)
    

    #multi_graph.GetXaxis().SetLimits(y_low_limit, y_max_limit)
    # Customize graph styles
    graph1_data.SetMarkerStyle(20)
    graph1_data.SetMarkerColor(ROOT.kBlack)
    graph1_data.SetLineColor(ROOT.kBlack)
    #graph1_data.GetXaxis().SetRangeUser(y_low_limit, y_max_limit)
    graph1_data.GetXaxis().SetLimits(y_low_limit, y_max_limit)
    graph2_mc.GetXaxis().SetLimits(y_low_limit, y_max_limit)

    
    graph2_mc.SetMarkerStyle(21)
    graph2_mc.SetMarkerColor(ROOT.kRed)
    graph2_mc.SetLineColor(ROOT.kRed)
    #graph2_mc.GetXaxis().SetRangeUser(y_low_limit, y_max_limit)
    # Create a TMultiGraph
    multi_graph = ROOT.TMultiGraph()
    # Add graphs to the TMultiGraph
    multi_graph.Add(graph1_data, "AP")  # "P" means draw points
    multi_graph.Add(graph2_mc, "AP")  # "L" means draw lines
    # Create a canvas
    canvas = ROOT.TCanvas("canvas", "Comparison", 900, 700)
    canvas.SetGrid()
    canvas.SetLogx()
    canvas.SetLeftMargin(0.24)
    #
    # Set axis titles
    multi_graph.SetTitle(f'{tittle};{x_axis_title};{y_axis_title}')
    # Draw the TMultiGraph
    multi_graph.Draw("AP")
    canvas.Update()
    


    # Adjust text sizes for axis titles and labels
    multi_graph.GetXaxis().SetTitleSize(0.04)  # Set X-axis title size
    multi_graph.GetYaxis().SetTitleSize(0.04)  # Set Y-axis title size
    multi_graph.GetXaxis().SetTitleOffset(1.2)
    multi_graph.GetXaxis().SetLimits(y_low_limit, y_max_limit)  # Fix X-axis range
    canvas.Modified()
    canvas.Update()
    #multi_graph.GetXaxis().SetLabelOffset(0.2)
    #
    
    #multi_graph.GetXaxis().SetRangeUser(0.0,1000)
    # Add a legend
    legend = ROOT.TLegend(0.9,0.6,1.0,0.8)
    legend.AddEntry(graph1_data, "DATA ", "pl")
    legend.AddEntry(graph2_mc, "MC", "pl")
    legend.SetFillStyle(0)  # 0 means no fill
    legend.SetBorderSize(0)  # Remove the border
    legend.Draw()
    
    # Save the canvas
    canvas.SaveAs(f"{file_name}.png")
   

#pt DGL
plot_comparison("Cosmics_muons_DATA_DGL.root", "Cosmics_muons_MC_DGL.root","ptmean_total", "ptmean_total", "", "p_{T} $\mu_{ref}$ [GeV]  ", "Mean of q/p_{T} relative residual  ",10, 1000.0, "DGL_Average_Resolution")
plot_comparison("Cosmics_muons_DATA_DGL.root", "Cosmics_muons_MC_DGL.root","ptSigma_total", "ptSigma_total", "", "p_{T} $\mu_{ref}$ [GeV]  ", "Width of q/p_{T} relative residual  ",10, 1000, "DGL_Sigma")

#dz DGL

#plot_comparison("Cosmics_muons_DATA_DGL.root", "Cosmics_muons_MC_DGL.root","dzmean_total", "dzmean_total", "", " |dz| ", "Mean of q/p_{T} relative residual  ",0.8, 200.0, "DGL_Average_Resolution_dz")
#plot_comparison("Cosmics_muons_DATA_DGL.root", "Cosmics_muons_MC_DGL.root","dzSigma_total", "dzSigma_total", "", " |dz| ", "Width of q/p_{T} relative residual  ",0.8, 200.0, "DGL_Sigma_dz")

#dxy DGL
#plot_comparison("Cosmics_muons_DATA_DGL.root", "Cosmics_muons_MC_DGL.root","dxymean_total", "dxymean_total", "", " |dxy| ", "Mean of q/p_{T} relative residual  ",0.8, 100.0, "DGL_Average_Resolution_dxy")
#plot_comparison("Cosmics_muons_DATA_DGL.root", "Cosmics_muons_MC_DGL.root","dxySigma_total", "dxySigma_total", "", " |dxy| ", "Width of q/p_{T} relative residual  ",0.8 , 100.0, "DGL_Sigma_dxy")


#pt DSA
plot_comparison("Cosmics_muons_DATA_DSA.root", "Cosmics_muons_MC_DSA.root","ptmean_total", "ptmean_total", "", "p_{T} $\mu_{ref}$ [GeV]  ", "Mean of q/p_{T} relative residual  ",10, 1000.0, "DSA_Average_Resolution")
plot_comparison("Cosmics_muons_DATA_DSA.root", "Cosmics_muons_MC_DSA.root","ptSigma_total", "ptSigma_total", "", "p_{T} $\mu_{ref}$ [GeV]  ", "Width of q/p_{T} relative residual  ",10, 1000, "DSA_Sigma")

#dz DSA

#plot_comparison("Cosmics_muons_DATA_DSA.root", "Cosmics_muons_MC_DSA.root","dzmean_total", "dzmean_total", "", " |dz| ", "Mean of q/p_{T} relative residual  ",0.8, 200.0, "DSA_Average_Resolution_dz")
#plot_comparison("Cosmics_muons_DATA_DSA.root", "Cosmics_muons_MC_DSA.root","dzSigma_total", "dzSigma_total", "", " |dz| ", "Width of q/p_{T} relative residual  ",0.8, 200.0, "DSA_Sigma_dz")

#dxy DSA
#plot_comparison("Cosmics_muons_DATA_DSA.root", "Cosmics_muons_MC_DSA.root","dxymean_total", "dxymean_total", "", " |dxy| ", "Mean of q/p_{T} relative residual  ",0.8, 100.0, "DSA_Average_Resolution_dxy")
#plot_comparison("Cosmics_muons_DATA_DSA.root", "Cosmics_muons_MC_DSA.root","dxySigma_total", "dxySigma_total", "", " |dxy| ", "Width of q/p_{T} relative residual  ",0.8 , 100.0, "DSA_Sigma_dxy")

