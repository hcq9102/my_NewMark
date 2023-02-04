#include <iostream>
#include <fstream>
#include <cmath>
#include <blaze/Blaze.h>
#include <string>
#define ToString(x)  (#x)

// how to convert variable name to string???????

//std::cout << "WR(2,.) data:  \n";
void getres(blaze::DynamicMatrix <double> v){
    std::string resname = ToString(v);
    std::ofstream fout("GS_func"+ resname+"_plot_res.csv");
    fout << "dt,WR\n";
    for (std::size_t step = 0; step <= 8000; step++){
        fout << (0.01)*step<< ","
             <<float(v(2,step)) << "\n";
    }
    fout.close();
}

// reimplement here:

// std::ofstream fout0("GS_WR_plot_res.csv");
//     fout0 << "dt,WR\n";
//     for (std::size_t step = 0; step <= nt; step++){
//         fout0 << dt*step<< ","
//              <<WR(2,step) << "\n";
//     }
//     fout0.close();

    // std::ofstream fout("Gaus_U_d_plot_res.csv");

    // fout << "dt,U_d\n";
    // for (std::size_t step = 0; step <= nt; step++){
    //     fout << dt*step<< ","
    //          <<U_d(2,step) << "\n";
    // }
    // fout.close();

    // std::ofstream fout1("Gaus_U_v_plot_res.csv");
    // fout1 << "dt,U_v\n";
    // for (std::size_t step = 0; step <= nt; step++){
    //     fout1 << dt*step<< ","
    //          <<U_v(2,step) << "\n";
    // }
    // fout1.close();

    // std::ofstream fout2("Gaus_U_a_plot_res.csv");
    // fout2 << "dt,U_a\n";
    // for (std::size_t step = 0; step <= nt; step++){
    //     fout2 << dt*step<< ","
    //          <<U_a(2,step) << "\n";
    // }
    // fout2.close();

    // // WR,WRa,WRv
    // // auto rs = rows( WR, { 2UL } );
    // // std::cout << "WR(2,.) data:  \n";
    // std::ofstream fout3("Gaus_WR_plot_res.csv");
    // fout3<< "dt,WR\n";
    // for (std::size_t step = 0; step <= nt; step++){
    //     fout3 << dt*step<< ","
    //         <<WR(2,step) << "\n";
    // }
    // fout3.close();

    // std::ofstream fout4("Gaus_WRv_plot_res.csv");
    // fout4 << "dt,WRv\n";
    // for (std::size_t step = 0; step <= nt; step++){
    //     fout4 << dt*step<< ","
    //          <<WRv(2,step) << "\n";
    // }
    // fout4.close();

    // std::ofstream fout5("Gaus_WRa_plot_res.csv");
    // fout5 << "dt,WRa\n";
    // for (std::size_t step = 0; step <= nt; step++){
    //     fout5<< dt*step<< ","
    //          <<WRa(2,step) << "\n";
    // }
    // fout5.close();
    

