# this is not nice... but it should work..

obsmap = [["Angularity"],["theta_g"]]
zcuts  = ["0.05","0.1","0.2"]
betas  = ["0","1","2"]
alphas = ["0.5","1","2"]

file=open("MC_LHJSS_Zjet.plot","w")

for  zcut in zcuts :
   for  beta in betas :
      for  alpha in alphas :
        file.write("# BEGIN PLOT /MC_LHJSS_Zjet/Angularity_zcut_"+zcut+"_beta_"+beta+"_alpha_"+alpha+"\n")
        file.write("Title=Angularity, $\\alpha$ = "+alpha+" $z_{cut}$ = "+zcut+" $\\beta$ = "+beta+"\n")
        file.write("LogY=0\n")
        file.write("XLabel=$e_\\alpha$\n")
        file.write("YLabel=$1/\\sigma \\text{d}\\sigma/\\text{d}e_{\\alpha}$\n")
        file.write("# END PLOT\n\n")

for  zcut in zcuts :
   for  beta in betas :
      for  alpha in alphas :

        file.write("# BEGIN PLOT /MC_LHJSS_Zjet/log_Angularity_zcut_"+zcut+"_beta_"+beta+"_alpha_"+alpha+"\n")
        file.write("Title=Angularity, $\\alpha$ = "+alpha+" $z_{cut}$ = "+zcut+" $\\beta$ = "+beta+"\n")
        file.write("LogX=1\n")
        file.write("LogY=0\n")
        file.write("XLabel=$e_\\alpha$\n")
        file.write("YLabel=$1/\\sigma \\text{d}\\sigma/\\text{d}e_\\alpha$\n")
        file.write("# END PLOT\n\n")


for  zcut in zcuts :
   for  beta in betas :
      for  alpha in alphas :

        file.write("# BEGIN PLOT /MC_LHJSS_Zjet/theta_g_zcut_"+zcut+"_beta_"+beta+"_alpha_"+alpha+"\n")
        file.write("Title=$\\theta_g$, $\\alpha$ = "+alpha+" $z_{cut}$ = "+zcut+" $\\beta$ = "+beta+"\n")
        file.write("LogY=0\n")
        file.write("XLabel=$\\theta_g$\n")
        file.write("YLabel=$1/\\sigma \\text{d}\\sigma/\\text{d}\\theta_g$\n")
        file.write("# END PLOT\n\n")

for  zcut in zcuts :
   for  beta in betas :
      for  alpha in alphas :

        file.write("# BEGIN PLOT /MC_LHJSS_Zjet/log_theta_g_zcut_"+zcut+"_beta_"+beta+"_alpha_"+alpha+"\n")
        file.write("Title=$\\theta_g$, $\\alpha$ = "+alpha+" $z_{cut}$ = "+zcut+" $\\beta$ = "+beta+"\n")
        file.write("LogX=1\n")
        file.write("LogY=0\n")
        file.write("XLabel=$\\theta_g$\n")
        file.write("YLabel=$1/\\sigma \\text{d}\\sigma/\\text{d}\\theta_g$\n")
        file.write("# END PLOT\n\n")






