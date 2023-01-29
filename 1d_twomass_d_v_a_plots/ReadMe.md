
1D_twomass_undamped folder:

Newmark data from U_dN,U_vN,U_aN from predictor_corrector:

        displacement res = U_dN
        
        velocity res = U_vN
        
        acceleration = U_aN
        
 Exact analytical solution:
 
        Appendix A
        The exact analytical solution to (9) can be expressed as:
        d1(t) =
        d2(t) = 
        
        implement the equations on the appendix A 方程组 d1(t) & d2(t);
        
        d1(t) will give you the displacement at time t for the first spring ?
        
        My question: how can I get the velocity, acceleration.
                     
               DO      differentiate d on t---> v
                       differentiate v on t---> a  
                     

NOTE:

        d,v,a plots source file is : 1D_twomass_undamp_dva_plot.ipynb  (saved on google drive as well)
        compare: Gauss seidel / Jacobi / Newmark / exact analytical solution
        
 
## Final Note: 
        final figures: 
        https://github.com/hcq9102/my_NewMark/tree/main/1d_twomass_d_v_a_plots/final%20plots%20of%20figure4
      
