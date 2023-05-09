/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "calibration.h"
#include "matrix_algo.h"


using namespace easy3d;



/**
 * TODO: Finish this function for calibrating a camera from the corresponding 3D-2D point pairs.
 *       You may define a few functions for some sub-tasks.
 * @return True on success, otherwise false. On success, the camera parameters are returned by
 */
bool Calibration::calibration(
        const std::vector<Vector3D>& points_3d, /// input: An array of 3D points.
        const std::vector<Vector2D>& points_2d, /// input: An array of 2D image points.
        double& fx,    /// output: focal length (i.e., K[0][0], which is equal to 'alpha' in our slides).
        double& fy,    /// output: focal length (i.e., K[1][1], which is equal to 'beta/sin(theta)' in our slides).
        double& cx,    /// output: x component of the principal point (i.e., K[0][2], which is 'u0' in our slides).
        double& cy,    /// output: y component of the principal point (i.e., K[1][2], which is 'v0' in our slides).
        double& skew,  /// output: skew factor (i.e., K[0][1], which is equal to '-alpha * cot(theta)' in our slides).
        Matrix33& R,   /// output: the 3x3 rotation matrix encoding camera rotation.
        Vector3D& t)   /// output：a 3D vector encoding camera translation.
{
    std::cout << "\nTODO: I am going to implement the calibration() function in the following file:\n"
                 "\t    - calibration_method.cpp\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tCamera calibration requires computing the SVD and inverse of matrices.\n"
                 "\tIn this assignment, I provide you with a 'Matrix' and a 'Vector' data structures for storing and\n"
                 "\tmanipulating matrices and vectors of arbitrary sizes. I also wrote some code to show you how to:\n"
                 "\t    - compute the SVD of a matrix;\n"
                 "\t    - compute the inverse of a matrix;\n"
                 "\t    - compute the transpose of a matrix.\n\n"
                 "\tFeel free to use any of the provided data structures and functions. The commonly used linear algebra\n"
                 "\tfunctions are provided in the following files:\n"
                 "\t    - Calibration/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Calibration/vector.h  Vectors of arbitrary dimensions and related functions.\n"
                 "\t    - Calibration/matrix_algo.h  Determinant, inverse, SVD, linear least-squares...\n"
                 "\tPlease refer to the above files for a complete list of useful functions and their usage.\n\n"
                 "\tIn your final submission, please\n"
                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
                 "\t    - include all the source code (and please do NOT modify the structure of the directories).\n"
                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
                 "\t    - make sure your code compiles and can reproduce your results without ANY modification.\n\n" << std::flush;


    std::cout << "\n[Liangliang]:\n"
                 "\tThe input parameters of this function are:\n"
                 "\t\t- points_3d: An array of 3D points (input to this function)\n"
                 "\t\t- points_2d: An array of 2D image points (input to this function)\n"
                 "\tThis function must return either 'true' on success or 'false' otherwise. On success, the camera\n"
                 "\tparameters are returned by the following variables:\n"
                 "\t\t- fx and fy: the focal lengths (in our slides, we use 'alpha' and 'beta')\n"
                 "\t\t- cx and cy: the principal point (in our slides, we use 'u0' and 'v0')\n"
                 "\t\t- skew:      the skew factor ('-alpha * cot_theta')\n"
                 "\t\t- R:         the 3x3 rotation matrix encoding camera orientation\n"
                 "\t\t- t:         a 3D vector encoding camera location.\n"
                 "\tIMPORTANT: don't forget to write your recovered parameters to the above variables." << std::endl;

    // TODO: check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)
    bool match = points_3d.size() == points_2d.size();
    if(match){
        if (points_2d.size() < 6) {
            std::cout << "incorrect size" << std::endl;
        }
    }
    else {
        std::cout << "incorrect size" << std::endl;
        return false;
    }

    // TODO: construct the P matrix (so P * m = 0).
    Matrix P(int(points_3d.size()) * 2, 12, 0.0);
    for (int i = 0; i < points_3d.size(); i++) {
        double u = points_2d[i].x();
        double v = points_2d[i].y();
        // set corresponding rows in P
        P.set_row(2 * i, { points_3d[i].x(),points_3d[i].y(),points_3d[i].z() ,1,0,0,0,0,
            -u * points_3d[i].x(),-u * points_3d[i].y(),-u * points_3d[i].z(),-u });
        P.set_row(2 * i+1, {0,0,0,0,points_3d[i].x(),points_3d[i].y(),points_3d[i].z() ,1,
            -v * points_3d[i].x(),-v * points_3d[i].y(),-v * points_3d[i].z(),-v });
    }
    // remember Pm = 0 and m is the flattened matrix of M = K[R|t]
    // TODO: solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
    int m = int(points_3d.size() * 2);
    int n = 12;
    Matrix U(m, m, 0.0);   // initialized with 0s
    Matrix S(m, n, 0.0);   // initialized with 0s
    Matrix V(n, n, 0.0); 
    // do svd decomposition and get M from V which subject to Pm = 0 and ||m|| = 1
    svd_decompose(P, U, S, V);
    
    Vector M =V.get_column(11);
    std::cout<<"SVD RESULT " << mult(P, M) <<"finished" << std::endl; // check if Pm = 0
    //   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
    //             should be very close to your input images points.

    // TODO: extract intrinsic parameters from M.
    Vector3D a3 = { M[8],M[9],M[10]};
    Vector3D a1 = { M[0],M[1],M[2]};
    Vector3D a2 = { M[4],M[5],M[6]};
    // calculate the value of rho
    double scale = 1 / a3.length();
    // output M and scale and a3
    std::cout << "scale: " << scale << std::endl;
    
    Vector a1xa3 = cross(a1, a3);
    Vector a2xa3 = cross(a2, a3);
    double Cosin_theta = - dot(a1xa3, a2xa3) / (a1xa3.length() * a2xa3.length());

    double sin_theta = sqrt(1 - Cosin_theta * Cosin_theta);
    
    double alpha = pow(scale,2) * a1xa3.length() * sin_theta;
    double beta = pow(scale,2) * a2xa3.length() * sin_theta;
    cx = pow(scale, 2) * dot(a1, a3);
    cy = pow(scale, 2) * dot(a2, a3);
    skew = -alpha * Cosin_theta / sin_theta;
    Matrix33 K(alpha, skew, cx,0,beta/sin_theta,cy,0,0,1);
    // K^-1
    Matrix invK;
    invK = inverse(K);
    Vector r1 =  cross(a2,a3)/a2xa3.length();
    Vector r3 = scale* a3;
    Vector r2 = cross(r3, r1);
    Vector3D b(M[3], M[7], M[11]);


    // extrinsic parameters
    t = scale * mult(invK, b);
    R = Matrix33(r1[0], r1[1], r1[2], r2[0], r2[1], r2[2], r3[0], r3[1], r3[2]);

    // TODO: extract extrinsic parameters from M.
   // test whether the corresponed point is correct or not, given the points_3d wether we obtain the points_2d
    Matrix extrinsic = Matrix(3, 4, 0.0); // should be like [R | t], which is extrinsic matrix
    extrinsic.set_column(0, R.get_column(0));
    extrinsic.set_column(1, R.get_column(1));
    extrinsic.set_column(2, R.get_column(2));
    extrinsic.set_column(3, t);
    // given a random 3d point from points_3d, we can obtain the 2d point, the result should be the same as the points_2d
    Vector3D test = points_3d[0];
    Vector4D test1 = { test[0],test[1],test[2],1 };
    // CALCULATE [SU,SV,S] = [K[R|t]] * [X,Y,Z,1], NOTICE THE TRUE 2D COORD SHOULD BE [SU/S,SV/S,S]
    Vector3D test2d = mult(K, mult(extrinsic, test1));
    // output the m and
    Matrix34 M1 = mult(K, extrinsic);
    // output the m and m1, m1 should be the same as SCLAE* m
    std::cout << "M1: " << M1 << std::endl;
    std::cout << "M: " << scale* M << std::endl;
    // out put and compare the 2d points ground truth and the points we obtain from the 3d points
    test2d /= abs(test2d[2]);
    std::cout << "test2d" <<" " << test2d << std::endl;
    

    std::cout << "points_2d[0]s" <<" " << points_2d[0] << std::endl;
    // decide the sign of scale
    if (test2d[2]-0 <-0.0001) {
        std::cout << "scale changed is: " << scale << std::endl;
        scale = -scale;

       
        
        std::cout << "scale changed is: " << scale << std::endl;
        // extrinsic parameters
        t = -t;
        R[2][0] *= -1;
        R[2][1] *= -1;
        R[2][2] *= -1;
    }
    fx = alpha;
    fy = beta/sin_theta;
    //    // Error analysis


    std::vector<double> errors_x;
    std::vector<double> errors_y;

    std::cout << std::left << std::setw(10) << "Point"
        << std::left << std::setw(20) << "Real x"
        << std::left << std::setw(20) << "Predicted x"
        << std::left << std::setw(15) << "Error x"
        << std::left << std::setw(20) << "Real y"
        << std::left << std::setw(20) << "Predicted y"
        << std::left << std::setw(15) << "Error y" << std::endl;

    std::cout << std::left << std::setw(10) << "----------"
        << std::left << std::setw(20) << "--------------------"
        << std::left << std::setw(20) << "--------------------"
        << std::left << std::setw(15) << "---------------"
        << std::left << std::setw(20) << "--------------------"
        << std::left << std::setw(20) << "--------------------"
        << std::left << std::setw(15) << "---------------" << std::endl;

    for (int i = 0; i < points_2d.size(); i++) {
        Vector3D temp3d = points_3d[i];
        Vector4D temp1 = { temp3d[0],temp3d[1],temp3d[2],1 };
        Vector3D temp2d = mult(K, mult(extrinsic, temp1));
        temp2d = temp2d / temp2d[2];

        double error_x = std::abs((points_2d[i][0] - temp2d[0]));
        double error_y = std::abs((points_2d[i][1] - temp2d[1]));

        errors_x.push_back(error_x);
        errors_y.push_back(error_y);

        std::cout << std::left << std::setw(10) << i
            << std::left << std::setw(20) << points_2d[i][0]
            << std::left << std::setw(20) << temp2d[0]
            << std::left << std::setw(15) << error_x
            << std::left << std::setw(20) << points_2d[i][1]
            << std::left << std::setw(20) << temp2d[1]
            << std::left << std::setw(15) << error_y << std::endl;
    }

    // calculate indicators
    double mse_x = 0.0;
    double mse_y = 0.0;
    double mae_x = 0.0;
    double mae_y = 0.0;
    double std_x = 0.0;
    double std_y = 0.0;

    for (int i = 0; i < errors_x.size(); i++) {
        mse_x += errors_x[i] * errors_x[i];
        mse_y += errors_y[i] * errors_y[i];
        mae_x += errors_x[i];
        mae_y += errors_y[i];
    }

    mse_x /= errors_x.size();
    mse_y /= errors_y.size();
    mae_x /= errors_x.size();
    mae_y /= errors_y.size();

    for (int i = 0; i < errors_x.size(); i++) {
        std_x += (errors_x[i] - mae_x) * (errors_x[i] - mae_x);
        std_y += (errors_y[i] - mae_y) * (errors_y[i] - mae_y);
    }

    std_x = std::sqrt(std_x / errors_x.size());
    std_y = std::sqrt(std_y / errors_y.size());

    // output error indicators
    std::cout << "MSE_x: " << mse_x << std::endl;
    std::cout << "MSE_y: " << mse_y << std::endl;
    std::cout << "MAE_x: " << mae_x << std::endl;
    std::cout << "MAE_y: " << mae_y << std::endl;
    std::cout << "STD_x: " << std_x << std::endl;
    std::cout << "STD_y: " << std_y << std::endl;



    
    std::cout << "\n\tTODO: After you implement this function, please return 'true' - this will trigger the viewer to\n"
                 "\t\tupdate the rendering using your recovered camera parameters. This can help you to visually check\n"
                 "\t\tif your calibration is successful or not.\n\n" << std::flush;
    return true;
}

















