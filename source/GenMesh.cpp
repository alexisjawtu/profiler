#include "GenMesh.h"


GenMesh::GenMesh(
    // to_delete const int type_,
    // to_delete const Point &pA_,
    // to_delete const Point &pD_,
    // to_delete double H_,
    // to_delete double h_,
    // to_delete double rescaled_diameter,
    // to_delete double theta_,
    // to_delete double L1_,
    // to_delete double L2_,
    // to_delete double angle_difference_,
    int n_layers,
    double h_layer,
    bool debug_flag,
    string prefix,
    map<string, valarray<double>>& profile_params
):
    
    // to_delete wtype(type_),
    // to_delete pA(pA_),
    // to_delete pD(pD_),
    // to_delete h_max(H_),
    // to_delete h_min(h_),
    // to_delete R(rescaled_diameter),
    // to_delete probe_theta(theta_),
    // to_delete L1(L1_),
    // to_delete L2(L2_),
    // to_delete angle_difference(angle_difference_),
    number_of_layers(n_layers),
    height_of_layer(h_layer),
    debug(debug_flag),
    output_dir(prefix)

{


    // TODO: remove these three, as we can call the map<> profile_params directly
    delta_over_plane_xy = profile_params["Width"];  // {delta_xy_0, delta_xy_1, delta_xy_2, delta_xy_3, delta_xy_4};
    upper_z             = profile_params["Ceiling"];  // {upper_z_0, upper_z_1, upper_z_2, upper_z_3, upper_z_4};
    lower_z             = profile_params["Floor"];  // {lower_z_0, lower_z_1, lower_z_2, lower_z_3, lower_z_4};

    case_translator.emplace(make_pair(false, false), make_pair(1, 0));
    case_translator.emplace(make_pair(false, true),  make_pair(3, 1));
    case_translator.emplace(make_pair(true,  false), make_pair(2, 0));
    case_translator.emplace(make_pair(true,  true),  make_pair(0, 1));

    outer_diagonals_translator.emplace(make_pair(false, false), make_pair(false, true));
    outer_diagonals_translator.emplace(make_pair(false, true), make_pair(false, false));
    outer_diagonals_translator.emplace(make_pair(true, false), make_pair(false, true));
    outer_diagonals_translator.emplace(make_pair(true, true), make_pair(true, false));

}


GenMesh::~GenMesh() {}


/*
double GenMesh::distance(const Point &p1, const Point &p2)
{
    double temp;
    temp = (p1.m_x - p2.m_x) * (p1.m_x - p2.m_x) + (p1.m_y - p2.m_y) * (p1.m_y - p2.m_y);
    return sqrt(temp);
}


double GenMesh::potential(const Point &p)
{
    return -log(distance(pA, p)) + log(distance(pD, p));
}


Point GenMesh::GradPotential(const Point& p, bool is_unit)
{
    Point temp;
    double tempx, tempy;
    double tempx2, tempy2;
    temp =  (pA - p) / (distance(p, pA) * distance(p, pA)) - (pD - p) / (distance(p, pD) * distance(p, pD));
    temp = temp * -1.0;
    if (is_unit)  temp = temp/temp.norm();
    return temp;
}


Point GenMesh::get_next_inital_p(double value)
{
    Point temp;
	double d = distance(pA, pD);
    double x;
	x = d/2.0*(1-exp(value))/(1+exp(value));
    temp = pAD + direc_AD*x;
    return temp;
}


Point GenMesh::get_next_inital_adjoint_p(double value)
{
    Point temp;
	double d = distance(pA, pD);
    double x;
    if( fabs(value) < 1E-6) return Point(0,1000); //The contour line goes to infinity.    
	x = - d/2.0*(1+exp(value))/(exp(value)-1);
    temp = pAD + direc_AD*x;
    return temp;
}

*/


Point GenMesh::get_right_orthogonal(Point vector2D)
{
    // TODO: test once more
    return Point(vector2D.m_y, -vector2D.m_x)/vector2D.norm();
}



// Point GenMesh::get_rotation_point(const Point& p, double theta, double h, const int direct)
// {
//     Point grad_direction = Point(0, 0) - GradPotential(p);
//     Point t;
//     t.m_x = grad_direction.m_y * direct;
//     t.m_y = -grad_direction.m_x * direct;
//     Point rotation_vec;
//     rotation_vec.m_x = (cos(theta) * t.m_x + sin(theta) * t.m_y) * h;
//     rotation_vec.m_y = (-sin(theta) * t.m_x + cos(theta) * t.m_y) * h;
//     Point temp = rotation_vec + p;
//     return temp;
// }


// Point GenMesh::get_next_p_on_contour(const Point &p0, double h, double value, const int direct)
// {
//     double last_theta = 0;
//     double cur_theta = -0.01 * PI;
//     int count = 0;
//     double cur_f_value = potential(get_rotation_point(p0, cur_theta*direct, h, direct)) - value;
//     double last_f_value = potential(get_rotation_point(p0, last_theta*direct, h, direct)) - value;
//     Point test1 = get_rotation_point(p0, cur_theta*direct, h, direct);
//     Point test2 = get_rotation_point(p0, last_theta*direct, h, direct);
//     double new_theta;
//     while (count < 8 && fabs(cur_theta - last_theta) > 1E-5)
//     {
//         new_theta = cur_theta - (cur_f_value) * (cur_theta - last_theta) / (cur_f_value - last_f_value);
//         PRINT(new_theta);
//         PRINT(cur_f_value);
//         PRINT(last_f_value);
//         last_theta = cur_theta;
//         cur_theta = new_theta;
//         last_f_value = cur_f_value;
//         cur_f_value = potential(get_rotation_point(p0, cur_theta*direct, h, direct)) - value;
//         count++;
//     }
//     return get_rotation_point(p0, new_theta*direct, h, direct);
// }


// int GenMesh::boundary_adjust_type1(Point& temp, double h, bool is_around_probe) {
//     // TODO: change name is_around_probe to is_in_the_inner_layer
//     // printf("(%.4f %.4f): R=%.3f,h=%.3f ", pCenter.m_x, pCenter.m_y, R, h);

//     double d = distance(pCenter, temp);
//     bool neighbor_point_check_on_boundary = true;

//     if (d < R - 0.3 * h) {
//         // The points falling in this case ought to be non boundary mesh points
//         return 1;
//     }

//     else {
//         if (d < R + 0.7 * h) {
//             double scale = R/d;
            
//             if (is_around_probe && d < R) {
//                 pointlist.push_back(temp);
//                 neighbor_point_check_on_boundary = false;
//             }

//             temp = pCenter + (temp - pCenter) * scale;

//             if (neighbor_point_check_on_boundary) {
//                 for (int i = 0; i < bdr_pointlist.size(); i++) {
//                     if (distance(temp, bdr_pointlist[i]) < h/3.0) {
//                         return 0;
//                     }
//                 }
//             }

//             /*
//                 push_to_bdr_pointlist:

//                     If it is true, we add the point to the bdr_pointlist.
//                     Remember bdr_pointlist[] and pointlist[] are two different lists,
//                     that is why we push the Point temp here to the bdr_pointlist and
//                     also in GenMesh::create_outer_points_type1() we push the same
//                     Point temp to pointlist, after the call to boundary_adjust_type1()
//                     is terminated.
//             */
//             bool push_to_bdr_pointlist = true;
//             int index = 0;

//             while (push_to_bdr_pointlist && index < bdr_pointlist.size()) {
//                 /*
//                     Traverse the present boundary points. If temp is "already one of them",
//                     discard it and don't push it to the boundary list
//                 */
//                 if (distance(temp, bdr_pointlist[index]) < minimum_float) {
//                     push_to_bdr_pointlist = false;
//                 }
//                 index ++;
//             }

//             if (push_to_bdr_pointlist) {
//                 bdr_pointlist.push_back(temp);
//             }

//             return 1;
//         }            
//     }
//     return 0;
// }


// int GenMesh::boundary_adjust_type2(Point &temp, double h, bool is_around_probe) {
//     // printf("(%.4f %.4f): R=%.3f,h=%.3f ", pCenter.m_x, pCenter.m_y, R,h );
//     bool neighbor_point_check_on_boundary = false;  // added 2022/3/30 (delete varible "flag")
//     double scale, arg; 
//     double d = distance(pCenter, temp);
//     arg = atan2(temp.m_y - pCenter.m_y, temp.m_x - pCenter.m_x);
//     if (arg < 0) arg += 2 * PI;

//     if (distance(temp, Orifla_start) < h) {
//         temp = Orifla_start;
//         neighbor_point_check_on_boundary = true;
//     } else if (distance(temp, Orifla_end) < h) {
//         temp = Orifla_end;
//         neighbor_point_check_on_boundary = true;  // modified 2022/3/30
//     }

//     if (neighbor_point_check_on_boundary) {  // modified 2022/3/30
//         for (int i = 0; i < bdr_pointlist.size(); i++) {
//             if (distance(temp, bdr_pointlist[i]) < h/3.0) {
//                 return 0;
//             }
//         }
//         bdr_pointlist.push_back(temp);
//         return 1;
//     }

//     if (angle_orifla_start < arg && arg < angle_orifla_end) {
//         double weight = cos(half_orifla_angle)/cos(arg - angle_orifla_mid);

//         if (d < (R - 0.3 * h) * weight) {
//             return 1;
//         } else {
//             if (d < (R + 0.9 * h) * weight) {
//                 scale = weight * R/d; 
//                 temp = pCenter + (temp - pCenter) * scale;
//                 neighbor_point_check_on_boundary = true;  // modified 2022/3/30

//                 if (distance(temp, Orifla_start) < h) {
//                     temp = Orifla_start;
//                 } else if (distance(temp, Orifla_end) < h) {
//                     temp = Orifla_end;
//                 } else {
//                     temp = temp + (temp - Orifla_mid) * 0.05;
//                 }

//             }
//         }
//     } else {
//         if (d < R - 0.3 * h) {
//             return 1;
//         } else {
//             if (d < R + 0.7 * h) {
//                 scale = R/d;
//                 neighbor_point_check_on_boundary = true;  // modified 2022/3/30
//                 temp = pCenter + (temp - pCenter) * scale;
//                 if (is_around_probe && d < R) {
//                     pointlist.push_back(temp);
//                     neighbor_point_check_on_boundary = false;
//                 }// added 2022/3/30 */
                
//             }            
//         }
//     }

//     if(neighbor_point_check_on_boundary){ // modified 2022/3/30
//         for(int i=0; i < bdr_pointlist.size(); i++){
//             if (distance(temp, bdr_pointlist[i])<h/3.0){
//                 return 0;
//             }
//         }
//         bdr_pointlist.push_back(temp);
//         return 1;
//     }

//     return 0;
// }


// int GenMesh::boundary_adjust_type3(Point &temp, double h, bool is_around_probe) {

//     // printf("(%.4f %.4f): R=%.3f,h=%.3f ", pCenter.m_x, pCenter.m_y, R,h );  
//     //int flag = 0;
//     bool neighbor_point_check_on_boundary = false;  //added 2022/3/30 (delete varible "flag")
//     double scale,arg;
//     double d= distance( pCenter, temp );
//     double delta = angle_2nd_orifla_end - 2 * PI;
//     arg = atan2(temp.m_y - pCenter.m_y, temp.m_x - pCenter.m_x);
//     // arg in [0, 2*pi)
//     if(arg< 0) arg += 2*PI;  
//     // if angle_2nd_orifla is gleater than 2*pi then arg in [delta, 2*pi + delta)
//     if(0 <= delta && arg < delta) arg += 2 * PI; 

//     if(distance(temp, Orifla_start) < h){
//         temp=  Orifla_start;
//         neighbor_point_check_on_boundary = true; //added 2022/3/30
//     }else if(distance(temp, Orifla_end) < h){
//         temp=  Orifla_end;
//         neighbor_point_check_on_boundary = true; //added 2022/3/30
//     }

//     if(distance(temp, Orifla_2nd_start) < h){
//         temp =  Orifla_2nd_start;
//         neighbor_point_check_on_boundary = true; //added 2022/3/30
//     }else if(distance(temp, Orifla_2nd_end) < h){
//         temp =  Orifla_2nd_end;
//         neighbor_point_check_on_boundary = true; //added 2022/3/30
//     }

//     if(neighbor_point_check_on_boundary == true ){ //modified 2022/3/30
//         for(int i=0; i < bdr_pointlist.size(); i++){
//             if (distance(temp, bdr_pointlist[i])<h/3.0){
//                 return 0;
//             }
//         }
//         bdr_pointlist.push_back(temp);
//         return 1;
//     }


//     if(angle_orifla_start < arg && arg <  angle_orifla_end ){
//         double weight = cos(half_orifla_angle)/cos(arg - angle_orifla_mid);
//         if( d < (R - 0.3*h)*weight ){
//             return 1;
//         }else{
//             if( d < (R + 0.9*h)*weight){
//                 scale = weight * R/d;
//                 temp = pCenter + (temp - pCenter) * scale;              
//                 temp =  temp + (temp - Orifla_mid) * 0.05;
//                 neighbor_point_check_on_boundary = true; //added 2022/3/30
//             }            
//         }
//     }else if(angle_2nd_orifla_start < arg && arg <  angle_2nd_orifla_end ){
//         double weight = cos(half_2nd_orifla_angle)/cos(arg - angle_2nd_orifla_mid);
//         if( d < (R - 0.3*h)*abs(weight) ){
//             return 1;
//         }else{
//             if (d < (R + 0.9*h)*abs(weight)) {
//                 scale = weight * R/d;
//                 temp = pCenter + (temp - pCenter) * scale;              
//                 temp =  temp + (temp - Orifla_2nd_mid) * 0.05;
//                 neighbor_point_check_on_boundary = true; //added 2022/3/30
//             }            
//         }
//     }else{
//         if( d < R - 0.3*h ){
//             return 1;
//         }else{
//             if( d < R + 0.7*h){
//                 scale = R/d;
//                 neighbor_point_check_on_boundary = true; //added 2022/3/30
//                 temp = pCenter + ( temp - pCenter) * scale;
//                 if(is_around_probe && d < R){
//                     pointlist.push_back(temp);
//                     neighbor_point_check_on_boundary = false;
//                 }// added 2022/3/30 */
                
//             }            
//         }
//     }

//     if (neighbor_point_check_on_boundary) {
//         for(int i = 0; i < bdr_pointlist.size(); i++){
//             if (distance(temp, bdr_pointlist[i]) < h/3.0) {
//                 return 0;
//             }
//         }
//         bdr_pointlist.push_back(temp);
//         return 1;
//     }

//     return 0;
// }


// int GenMesh::boundary_adjust_type4(Point &temp, double h, bool is_around_probe){

//     // printf("(%.4f %.4f): R=%.3f,h=%.3f ", pCenter.m_x, pCenter.m_y, R,h );  
//     //int flag = 0;
//     bool neighbor_point_check_on_boundary = false;  //added 2022/3/30
//     double scale,arg,r0,notch_angle,dist_ep;
//     Point Notch_start, Notch_end;
//     notch_angle = 0.5 * (angle_notch_end - angle_notch_start);
//     r0 = R*cos(notch_angle) - Foffset;
//     Notch_start = Point(-0.5*Owidth, R*sin(angle_notch_start));
//     Notch_end = Point(0.5*Owidth, R*sin(angle_notch_end));

//     double d= distance( pCenter, temp );
//     arg = atan2(  temp.m_y - pCenter.m_y, temp.m_x - pCenter.m_x );
//     if(arg< 0) arg += 2*PI;

//     if(distance(temp, Notch_start) < h){
//         temp=  Notch_start;
//         neighbor_point_check_on_boundary = true;  //added 2022/3/30
//     }else if(distance(temp, Notch_end) < h){
//         temp=  Notch_end;
//         neighbor_point_check_on_boundary = true;  //added 2022/3/30
//     }

//     if(neighbor_point_check_on_boundary == true){ //modified 2022/3/30
//         for(int i=0; i < bdr_pointlist.size(); i++){
//             if (distance(temp, bdr_pointlist[i])<h/3.0){
//                 return 0;
//             }
//         }
//         bdr_pointlist.push_back(temp);
//         return 1;
//     }

//     if(angle_notch_start < arg && arg < angle_notch_end){
//         if( d < R - L2 - 0.3*h ){ 
//             return 1;
//         }else{
//             return 0;
//         }
//     }else{
//         if( d < R - 0.3*h ){ 
//             return 1;
//         }else{
//             if( d < R + 0.7*h){
//                 scale = R/d;
//                 neighbor_point_check_on_boundary = true;  //added 2022/3/30
//                 temp = pCenter + ( temp - pCenter) * scale;
//                 if(is_around_probe && d < R){
//                     pointlist.push_back(temp);
//                     neighbor_point_check_on_boundary = false;
//                 }// added 2022/3/30 */
                
//             }            
//         }
//     }

//     if(neighbor_point_check_on_boundary == true){ //modified 2022/3/30
//         for(int i=0; i < bdr_pointlist.size(); i++){
//             if (distance(temp, bdr_pointlist[i])<h){
//                 return 0;
//             }
//         }
//         bdr_pointlist.push_back(temp);
//         return 1;
//     }
//     return 0;
// }


// bool GenMesh::is_inner_layer(Point &temp, double h) {

//     // printf("(%.4f %.4f): R=%.3f,h=%.3f ", pCenter.m_x, pCenter.m_y, R,h );     

//     double d = distance(pAD, temp);
//     int flag = 0;
//     if (d < innerR - 0.25 * h) {
//         return true;
//     } else {
//         if (d < innerR + 0.75 * h) {
//             double scale = innerR/d;
//             temp = pAD + (temp - pAD) * scale;
//             //return true;
//             flag = 1;
//         }
//     }

//     if (flag == 1) {
//         for (Point p: pointlist) {
//             /*
//             Check if it already belongs to pointlist
//             */
//             if (distance(temp, p) < 1e-8) {
//                 return false;
//             }
//         }
//         return true;
//     }
//     return false;
// }


// void GenMesh::get_dist_probe_to_bound_type1() {

//     dist_probe_to_bound_max = 0.0;
//     dist_probe_to_bound_min = 2.0*R;
//     int num = 100;
//     double dt = 2.0*PI/num, arg, d;
//     Point temp;

//     for (int i = 0; i < num; i++) {
//         arg = i * dt;
//         temp = pCenter + Point(cos(arg), sin(arg)) * R;
//         d = distance(temp, pAD);
//         if (d > dist_probe_to_bound_max) dist_probe_to_bound_max = d;
//         if (d < dist_probe_to_bound_min) dist_probe_to_bound_min = d;
//     }

// }


// void GenMesh::get_dist_probe_to_bound_type2(){
//     dist_probe_to_bound_max = 0.0;
//     dist_probe_to_bound_min = 2.0*R;

//     int num = 100;
//     double dt = 2.0 * PI/num, arg, d;
//     Point temp;
//     for (int i = 0; i < num; i++) {
//         arg = i * dt;
//         if (angle_orifla_start < arg && arg < angle_orifla_end) {
//             temp = pCenter + Point(cos(arg), sin(arg)) * R * cos(half_orifla_angle)/cos(arg - angle_orifla_mid);
//             d = distance(temp, pAD);
//         } else {
//             temp = pCenter + Point(cos(arg), sin(arg)) * R;
//             d = distance(temp, pAD);
//         }
//         if (d > dist_probe_to_bound_max) dist_probe_to_bound_max = d;
//         if (d < dist_probe_to_bound_min) dist_probe_to_bound_min = d;
//     }

// }


// void GenMesh::get_dist_probe_to_bound_type3(){
//     dist_probe_to_bound_max = 0.0;
//     dist_probe_to_bound_min = 2.0*R;
//     int num = 100;
//     double dt = 2.0*PI/num, arg, delta, d;
//     Point temp;
//     delta = angle_2nd_orifla_end - 2 * PI;
//     for(int i=0; i<num; i++){
//         arg = i * dt;
//          // if angle_2nd_orifla is gleater than 2*pi then arg in [delta, 2*pi + delta)
//         if(0 <= delta && arg < delta) arg += 2 * PI; 
//         if(angle_orifla_start < arg && arg < angle_orifla_end){
//             temp = pCenter + Point(cos(arg), sin(arg)) * R * cos(half_orifla_angle)/ cos(arg - angle_orifla_mid);
//             d = distance(temp, pAD);
//         }else if(angle_2nd_orifla_start < arg && arg < angle_2nd_orifla_end){
//             temp = pCenter + Point(cos(arg), sin(arg)) * R * cos(half_2nd_orifla_angle)/ cos(arg - angle_2nd_orifla_mid);
//             d = distance(temp, pAD);
//         }else{
//             temp = pCenter + Point(cos(arg), sin(arg))*R;
//             d = distance(temp, pAD);
//         }
//         if(d>dist_probe_to_bound_max)  dist_probe_to_bound_max = d;
//         if(d<dist_probe_to_bound_min)  dist_probe_to_bound_min = d;
//     }

// }


// void GenMesh::get_dist_probe_to_bound_type4(){
//     dist_probe_to_bound_max = 0.0;
//     dist_probe_to_bound_min = 2.0*R;

//     int num = 200;
//     double dt = 2.0*PI/num, arg, d, dist_ep, notch_angle,r0;
//     Point temp, notch_start, notch_end;
//     notch_angle = 0.5 * (angle_notch_end - angle_notch_start);
//     r0 = R*cos(notch_angle) - Foffset;
//     notch_start = Point(-0.5*Owidth, R*sin(angle_notch_start));
//     notch_end = Point(0.5*Owidth, R*sin(angle_notch_end));

//     for(int i=0; i<num; i++){
//         arg = i * dt;
//         if(angle_notch_start < arg && arg < angle_notch_end){            
//             temp = pCenter + Point(cos(notch_angle)/sin(arg)*cos(arg), sin(angle_notch_start)) * R ;
//             if(-L2/sqrt(2) <  temp.m_x && temp.m_x < L2/sqrt(2) ){
//                 dist_ep = abs(temp.m_x - (notch_start.m_x + notch_end.m_x)*0.5 );
//                 temp.m_y = sqrt(L2*L2-dist_ep*dist_ep) - r0;
//             }else if(temp.m_x < 0 ){
//                 dist_ep = abs(temp.m_x - notch_start.m_x);
//                 temp.m_y = temp.m_y + dist_ep;
//             }else{             
//                 dist_ep = abs(temp.m_x - notch_end.m_x);
//                 temp.m_y = temp.m_y + dist_ep;
//             }
//             d = distance(temp, pAD);

//         }else{   
//             temp = pCenter + Point( cos(arg), sin(arg))*R;
//             d = distance(temp, pAD);
//         }
//         if (d > dist_probe_to_bound_max) dist_probe_to_bound_max = d;
//         if (d < dist_probe_to_bound_min) dist_probe_to_bound_min = d;
//     }
 
// }


// void GenMesh::create_inner_points()
// {

//     double value; // value for contour lines
//     double step;  // each step for value variation.
//     int step_num;
//     h_min = h_min/4.0;
//     Point p = pA + direc_AD*2*h_min;
//     Point grad = GradPotential(p, false);

//     double base_value = potential(p);
//     double grad_norm = grad.norm();
//     double local_h = h_min;
//     double scale = local_h * grad_norm;
//     double last_value;
//     double value_B;

//     Point current;

//     double local_h_min=10000.0;
//     double r_min = 10000.0;
//     double r;
//     double value_max, value_min=1000;
//     double AD_distance = distance(pA, pD);
//     vector<double> value_list;

//     // First point of the mesh
//     pointlist.push_back(pA);
//     // Second point of the mesh   
//     pointlist.push_back(pD);
    
//     inner_max_h = h_min*2;

//     value = value_min;

//     //Calculate the step numbers. 
//     double _v1, _v2, _dv ;
//     p = pA + direc_AD*2*h_min; _v1 = potential(p);
//     p = pA + direc_AD*3*h_min; _v2 = potential(p);
//     step = _v1  - _v2;

//     p = pA + direc_AD*(2.0*h_min);
//     value_max = potential(p);
//     value_min = -potential(p);

//     p = pA + direc_AD*dist_AD/3.0;
//     value_B = potential(p);

//     std::cout << "Point A: (" << pA <<")\n";
//     std::cout << "Point B: (" << p  <<")\n";

//     step_num = round((value_max - value_B)/step+0.5);
//     step = (value_max - value_B)/step_num;
//     //Value of potential from A+2*h_min to B-epsilon.
//     value = value_max;
//     while(value > value_B + 1E-10) {
//         value_list.push_back(value);
//         value = value - step;
//     }

//     step_num = round(value_B/step+0.5);
//     step = value_B/step_num;
//     //Value of potential from B-epsilon to C+epsilon.
//     value = value_B;
//     while(value > -value_B - 1E-10){
//         value_list.push_back(value);
//         value = value - step;
//     }

//     //Value of potential from C+epsilon to D-2*h_min.
//     step_num = round((value_max - value_B)/step+0.5);
//     step = (value_max - value_B)/step_num;
//     value = -value_B - step;
//     while(value > value_min - 1E-10){
//         value_list.push_back(value);
//         value = value - step;
//     }

//     // Leave here for debug
//     // std::cout << "Min value and Max value, Base value, Scale " <<  value_min << " "  << value_max <<" "
//     //           << base_value <<" "  << scale <<"\n";
//     // std::cout << "Step around probes, Step number " <<  step << " "  << step_num <<" " <<"\n";

//     for (int i = 0; i < value_list.size(); i++) {
//         value = value_list[i];

//         Point initial = get_next_inital_p(value);
//         Point initial_adjoint = get_next_inital_adjoint_p(value);
//         double dist_to_adjoint;
        
//         r = distance(pA, initial);
//         if (r < r_min) {r_min = r;}
        
//         grad = GradPotential(initial, false);
//         grad_norm = grad.norm();
//         double h0 = scale / grad_norm;
//         local_h = h0;

//         if (is_inner_layer(initial, h0)) {
//             pointlist.push_back(initial);
//         }
               
// 		// Above part 
//         int direct_list[] = {-1, 1};
//         for(int direct_idx = 0; direct_idx < 2; direct_idx++) {

//             current = initial;
//             int count = 0;
//             dist_to_adjoint = distance(current, initial_adjoint);

//             while (dist_to_adjoint > local_h) //Start from initial point with direct_list[idx] for ONLY half cricle.
//             {
//                 grad = GradPotential(current, false);
//                 local_h = scale / grad.norm();
//                 if(local_h > h_max) local_h =h_max;                
//                 if(local_h < local_h_min) local_h_min =local_h;
                
//                 current = get_next_p_on_contour(current, local_h, value, direct_list[direct_idx]);
//                 dist_to_adjoint = distance(current, initial_adjoint);

//                 if (dist_to_adjoint < local_h/3.0) break;
                
//                 bool flag_is_inner_layer = is_inner_layer(current, local_h);
//                 if (flag_is_inner_layer) {
//                     // (2022/1/1)  "local_h*2.7" comes from that this line move to boundary or delete neighboring
//                     // point x satisfying dist(BoundaryOmega, x) < local_h * 0.9 (see, e.g. line 187). 
//                     if (wtype == 1 && boundary_adjust_type1(current, local_h*2.7, flag_is_inner_layer)) {
//                         pointlist.push_back(current);
//                     }

//                     if (wtype == 2 && boundary_adjust_type2(current, local_h * 2.7, flag_is_inner_layer)) {
//                         pointlist.push_back(current); 
//                     }

//                     if (wtype == 3 && boundary_adjust_type3(current, local_h * 2.7, flag_is_inner_layer)) {
//                         pointlist.push_back(current);
//                     }

//                     if (wtype == 4 && boundary_adjust_type4(current, local_h * 0.9, flag_is_inner_layer)) {
//                         pointlist.push_back(current);    
//                     }

//                 } else {
//                     break;
//                 }

//                 if (inner_max_h < local_h && direc_AD * (current - pA) < 0) inner_max_h = local_h;
//                 if (inner_max_h < local_h) inner_max_h = local_h;

//                 count ++;
//                 if (count > 1000) break;
//             }
        
//             if (direct_idx > 0) continue;
            
//             Point temp = get_next_inital_adjoint_p(value);
//             bool flag_is_inner_layer = is_inner_layer(temp, local_h);
            
//             if (flag_is_inner_layer) {
//                 if (wtype == 1 && boundary_adjust_type1(temp, local_h * 2.7, flag_is_inner_layer)) {
//                     pointlist.push_back(temp);
//                 }
//                 if (wtype == 2 && boundary_adjust_type2(temp, local_h * 2.7, flag_is_inner_layer)) {
//                     pointlist.push_back(temp);
//                 }
//                 if (wtype == 3 && boundary_adjust_type3(temp, local_h * 2.7, flag_is_inner_layer)) {
//                     pointlist.push_back(temp);
//                 }
//                 if (wtype == 4 && boundary_adjust_type4(temp, local_h * 0.9, flag_is_inner_layer)) {
//                     pointlist.push_back(temp);
//                 }
//             }
//         }
//     }

//     //Add outer circle points.
    
//     std::cout << "Value of variable inner_max_h: " << inner_max_h << "\n";

//     r_min = r_min *0.5;

//     int div_num=9*2;
//     double dt=2*PI/div_num;

//     for (int k = 0; k < div_num; k++) {
//         pointlist.push_back(pD + Point(cos(dt * k + probe_theta), sin(dt * k + probe_theta)) * r_min); 
//         pointlist.push_back(pA + Point(cos(dt * k + probe_theta), sin(dt * k + probe_theta)) * r_min);
//     }

//     cout << "Size of inner layer " << pointlist.size() << endl;

// }


// void GenMesh::create_outer_points_type1(){
//     double local_max_h = h_max;
//     double local_min_h = inner_max_h;// *.5; //2022/3/30 add *.5
//     double dist_from_center =0;
//     double h, theta;
    
//     std::cout <<"local max h: " << local_max_h  << " | local min h: " << local_min_h << "\n";

//     Point initial = Point(0, innerR + local_min_h);
//     Point temp;
//     int tmp_k=-1;
//     dist_from_center = innerR + local_min_h;
//     h = local_min_h;

//     while( dist_from_center < 2*R+h){
//         int num = 2.0*PI*dist_from_center / h;
//         double d_theta = 2.0*PI/num;
//         double t = ( dist_from_center - innerR )/ (dist_probe_to_bound_max - innerR);
//         h = local_max_h * t + local_min_h*(1-t);

//         for(int i=0; i<num; i++){
//             theta = i*d_theta + probe_theta;
//             if( tmp_k == 1) theta += d_theta/2;
//             temp = pAD + Point( dist_from_center*cos(theta), dist_from_center*sin(theta));
//             if(boundary_adjust_type1(temp,h)) {
//                 pointlist.push_back(temp);
//             }
//         }

//         dist_from_center += h; 
//         tmp_k *=-1;
//     }
//     std::cout << "probe_theta: " << probe_theta <<"\n";
// }


// void GenMesh::create_outer_points_type2() {

//     double local_max_h = h_max;
//     double local_min_h = inner_max_h;  // *.5; //2022/3/30 add *.5
//     double dist_from_center = 0;
//     double h, theta;

//     std::cout << "local max h: " << local_max_h  << " | local min h: " << local_min_h << "\n";

//     Point initial = Point(0, innerR + local_min_h);
//     Point temp;
//     int tmp_k = -1;
//     dist_from_center = innerR + local_min_h;
//     h = local_min_h;

//     pointlist.push_back(Point(cos(angle_orifla_start), 0));
//     while (dist_from_center < 2 * R + h) {
//         int num = 1 + 2.0 * PI * dist_from_center/h;
//         double d_theta = 2.0 * PI/num;
//         double t = (dist_from_center - innerR)/(dist_probe_to_bound_max - innerR);
//         h = local_max_h * t + local_min_h * (1 - t);

//         for (int i = 0; i < num; i++) {
//             theta = i * d_theta + probe_theta;
//             if (tmp_k == 1) theta += d_theta/2;
//             temp = pAD + Point(dist_from_center * cos(theta), dist_from_center * sin(theta));
//             if (boundary_adjust_type2(temp, h)) {
//                 pointlist.push_back(temp);
//             }
//         }

//         dist_from_center += h; 
//         tmp_k *= -1;
//     }

//     std::cout << "probe_theta: " << probe_theta <<"\n";

// }


// void GenMesh::create_outer_points_type3(){
//     double local_max_h = h_max;
//     double local_min_h = inner_max_h;//*.5; //2022/3/30 add *.5
//     double dist_from_center =0;
//     double h, theta;
    
    
//     std::cout <<"local max h: " << local_max_h  << " local min h " << local_min_h << "\n";

//     Point initial = Point(0, innerR + local_min_h);
//     Point temp;
//     int tmp_k=-1;
//     dist_from_center = innerR + local_min_h;
//     h = local_min_h;

//     while( dist_from_center < 2*R+h){
//         int num = 1 + 2.0*PI*dist_from_center / h;
//         double d_theta = 2.0*PI/num;
//         double t = ( dist_from_center - innerR )/ (dist_probe_to_bound_max - innerR);
//         h = local_max_h * t + local_min_h*(1-t);

//         for(int i=0; i<num; i++){
//             theta = i*d_theta + probe_theta;
//             if( tmp_k == 1) theta += d_theta/2;
//             temp = pAD + Point( dist_from_center*cos(theta), dist_from_center*sin(theta));

//             if(boundary_adjust_type3(temp,h)){
//                 pointlist.push_back( temp );
//             }
//         }

//         dist_from_center += h; 
//         tmp_k *=-1;
//     }
//     std::cout << "probe_theta " << probe_theta <<"\n"; 


// }


// void GenMesh::create_outer_points_type4() {
//     double dist_ep, notch_angle,r0;
//     double theta, dt;
//     Point temp, notch_start, notch_end;
//     notch_angle = 0.5 * (angle_notch_end - angle_notch_start);
//     r0 = R*cos(notch_angle) - Foffset;
//     notch_start = Point(-0.5*Owidth, R*sin(angle_notch_start));
//     notch_end = Point(0.5*Owidth, R*sin(angle_notch_end));

//     //initR = R;
//     //R = r0;   
//     create_outer_points_type1();

//     //Adjustment for all boundary points
//     //rmove notch part(Point list) 
//     for(auto p = pointlist.begin(); p != pointlist.end();){
//         temp=*p;
//         double arg = atan2(temp.m_y, temp.m_x) ;
//         if(arg < 0) arg += 2*PI;
//         if( angle_notch_start - 1e-3 <= arg && arg <= angle_notch_end + 1e-3 && abs(temp.norm()-R) < L1){
//             p = pointlist.erase(p);
//         }else{
//             ++p;
//         }
//     }

//     //rmove notch part(Boundary point list)
//     for(auto p = bdr_pointlist.begin(); p != bdr_pointlist.end();){
//         temp=*p;
//         double arg = atan2(temp.m_y, temp.m_x) ;
//         if(arg < 0) arg += 2*PI;
//         if( angle_notch_start - 1e-3 <= arg && arg <= angle_notch_end + 1e-3 && abs(temp.norm()-R) < L1){
//             p = bdr_pointlist.erase(p);
//         }else{
//             ++p;
//         }
//     }


//     //Notch shape
//     int num = 10;
//     dt = 2.0*notch_angle/num;
//     theta = angle_notch_start ;
//     vector<Point> notch_pointlist;
//     while(1){            
//         temp = pCenter + Point(cos(notch_angle)/sin(theta)*cos(theta), sin(angle_notch_start)) * R ;

//         if(abs(theta - angle_notch_start) < 1e-6){
//             bdr_pointlist.push_back( temp ); 
//             notch_pointlist.push_back( temp );
//          }
        
//         if(-L2/sqrt(2) <  temp.m_x && temp.m_x < L2/sqrt(2) ){
//             dist_ep = abs(temp.m_x - (notch_start.m_x + notch_end.m_x)*0.5 );
//             temp.m_y = sqrt(L2*L2-dist_ep*dist_ep) - r0;
//             bdr_pointlist.push_back( temp ); 
//             notch_pointlist.push_back( temp );
//         }

//         if(abs(theta - angle_notch_end) < 1e-6){
//             bdr_pointlist.push_back( temp ); 
//             notch_pointlist.push_back( temp );
//             break;
//         }
//         theta += dt;
//     } 
  
//     //Find nearest point from notch.
//     double s_radius = 2.0 * R ;
//     for(Point p: pointlist){
//         for(Point q: notch_pointlist){
//             if( 1e-8 < distance( p, q ) && 1e-3 < abs(p.norm()-R) && distance( p, q ) < s_radius  ){
//                 s_radius = distance( p, q );
                 
//             }
//         }
//     }
//     cout << "dist(inner points, notch): " << s_radius << endl;

//     //Connect notch part.
//     pointlist.insert(pointlist.end(), notch_pointlist.begin(), notch_pointlist.end()); 


//     //Around notch
//     s_radius *= 0.8;
//     int r_div_num = s_radius/L1;
//     int a_div_num = (num + 1) - r_div_num;
//     while(L1 < s_radius){
//         dt = PI/a_div_num; 
//         for(int i=0; i<=a_div_num; i++){
//             theta = i*dt;
//             if(a_div_num%2 == 0) theta += 0.25*PI;
//             temp = Point(cos(theta), sin(theta))*s_radius;
//             temp = temp + Point(0, -1.0*R);
            
//             if( boundary_adjust_type4(temp, L1) ){
//                 pointlist.push_back( temp );
//             } 
    
//         }
//         a_div_num += 1;
//         s_radius -= L1;
//     }
    
//     //Exterior point (for removing notch part)
//     pointlist.push_back( Point(0, -1.0*R) );

// }


// int GenMesh::homogenize_point_distribuition_around_probe(double radious){

//     Point bdrpoint_of_pCenter2pAD = pCenter + (pAD-pCenter) * (R/distance(pAD, pCenter));

//     if(wtype == 2){
//         if(angle_orifla_start < arg_pAD && arg_pAD < angle_orifla_end) return 0;
//     }else if(wtype==3){
//         if(angle_orifla_start < arg_pAD && arg_pAD < angle_orifla_end) return 0;
//         double tmp = acos( pAD*Orifla_2nd_mid/pAD.norm()/Orifla_2nd_mid.norm() ); // in [0, 0.5*pi]
//         if(tmp < fabs(angle_2nd_orifla_mid -angle_2nd_orifla_start)) return 0;
//     }else if(wtype==4){
//         if(angle_notch_start < arg_pAD && arg_pAD < angle_notch_end) return 0;
//     }

//     bool add_bdrpoint_of_pCenter2pAD = true;

//     for (auto idx = bdr_pointlist.begin(); idx != bdr_pointlist.end();) {
//         Point temp = *idx;
//         if(distance(temp, bdrpoint_of_pCenter2pAD) < 1e-8 ) { 
//             add_bdrpoint_of_pCenter2pAD = false; 
//             break;
//         }
//         ++idx;
//     }

//     if (add_bdrpoint_of_pCenter2pAD) {
//         /*

//             IMPORTANT FOR TYPE 4:
                
//                 pointlist.pushback doesn't work for type 4!!

//                 Replace the following line with:

//                     pointlist.insert(pointlist.end() - 1, bdrpoint_of_pCenter2pAD); 

//         */
//         pointlist.push_back(bdrpoint_of_pCenter2pAD);

//         bdr_pointlist.push_back(bdrpoint_of_pCenter2pAD);
//     }

//     return 0;

// }


// vector<double> GenMesh::convert_double_lst() {
//     vector<double> datalist;
    
//     for(size_t i = 0; i < pointlist.size(); ++i) {
//         datalist.push_back(pointlist[i].m_x);
//         datalist.push_back(pointlist[i].m_y);
//     }
    
//     return datalist;

// }


void GenMesh::stream_elements_out() {
    ofstream outfile;
    vector<int> current;
    outfile.open(output_dir + "/profile_elements.dat", ios::trunc);

    for (int e = 0; e < elements_by_vertices.size(); e++) {
        current = elements_by_vertices[e];
        outfile << current[0] << ", " << current[1] << ", "
                << current[2] << ", " << current[3] << endl;
    }
    outfile.close();
}


// TODO: use a template here.
void GenMesh::print_vector(vector<int>& v) {
    for (auto& d: v) {
        cout << d << ", ";
    }
    cout << "\n";
}


void GenMesh::print_vector(vector<std::size_t>& v) {
    for (auto& d: v) {
        cout << d << ", ";
    }
    cout << "\n";
}


void GenMesh::print_vector(vector<double>& v) {
    for (auto& d: v) {
        cout << d << ", ";
    }
    cout << "\n";
}


void GenMesh::print_vector(vector<bool>& v) {
    for (int i = 0; i < v.size(); i++) {
        cout << boolalpha << v[i] << ", ";
    }
    cout << "\n";
}


Point GenMesh::project2D(Point3D p) {
    return Point(p.x, p.y);
}


void GenMesh::radial_sort_boundary_points() {
    vector<Point> original_points(bdr_pointlist);

    int minimal_index;
    Point minimal_point;

    for (int j = 0; j < bdr_pointlist.size(); j++) {

        minimal_index = 0;
        minimal_point = original_points[0];

        for (int r = 1; r < original_points.size(); r++) {
            if (original_points[r] < minimal_point) {
                minimal_point = original_points[r];
                minimal_index = r;
            }
        }

        sorted_boundary_points.emplace(j, Point(minimal_point));
        radial_sorted_upper_bdr_indices.push_back(upper_bdr_points_global_indices[minimal_index]);

        original_points.erase(original_points.begin() + minimal_index);
        upper_bdr_points_global_indices.erase(upper_bdr_points_global_indices.begin() + minimal_index);
    }

    cuadrilaterals = radial_sorted_upper_bdr_indices.size();

    for (int j = 0; j < cuadrilaterals; j++)
        radial_sorted_lower_bdr_indices.push_back(pointlist.size() + radial_sorted_upper_bdr_indices[j]);

    // Repeat first indices to have the last cuadrilateral.
    if (!radial_sorted_upper_bdr_indices.size())
    {
        cout << "Trouble with filling radial_sorted_upper_bdr_indices.\n"
                "In function void GenMesh::radial_sort_boundary_points." << endl;
        exit(1);
    }
    radial_sorted_upper_bdr_indices.push_back(radial_sorted_upper_bdr_indices[0]);
    radial_sorted_lower_bdr_indices.push_back(radial_sorted_lower_bdr_indices[0]);

}



// void GenMesh::new_radial_sort_boundary_points() {
//     /**

//     *   This is on testing stage.

//     *   After that, we replace the former radial_sort_boundary_points()

//     */
//     vector<Point> original_points(bdr_pointlist);

//     int fixed_size = original_points.size();
//     int minimal_index;
//     Point minimal_point;

//     for (int j = 0; j < fixed_size; j++) {

//         minimal_index = 0;
//         minimal_point = original_points[0];

//         for (int r = 1; r < original_points.size(); r++) {
//             if (original_points[r] < minimal_point) {
//                 minimal_point = original_points[r];
//                 minimal_index = r;
//             }
//         }

//         sorted_boundary_points.emplace(j, Point(minimal_point));
//         radial_sorted_upper_bdr_indices.push_back(upper_bdr_points_global_indices[minimal_index]);

//         original_points.erase(original_points.begin() + minimal_index);
//         upper_bdr_points_global_indices.erase(upper_bdr_points_global_indices.begin() + minimal_index);
//     }
// }



// vector<Point> GenMesh::get_boundary_points() const {
//     return bdr_pointlist;
// }


void GenMesh::orient_profile_diagonals() {
    /*
    
    Cases of diagonals. 

    Viewing the profile walls from <<outside>> the wafer, those are:

                 o------o                  o------o
        case 0:  |   /  |         case 1:  |  \   |
                 |  /   |                  |   \  |
                 o------o                  o------o
    
    */

    if (debug) {
        cout << "Number of vertical cuadrilaterals around profile: " 
             << cuadrilaterals << endl;
    }

    for (auto& elements_list: vertices_by_elements)
        elements_list.second.sort();

    vector<int>::iterator a;
    list<int> first, second;
 
    for (int c = 0; c < cuadrilaterals; c++) {
        first  = vertices_by_elements[radial_sorted_upper_bdr_indices[c]];
        second = vertices_by_elements[radial_sorted_lower_bdr_indices[c + 1]];

        vector<int> intersection(max(first.size(), second.size()));
    
        a = set_intersection(first.begin(), first.end(), second.begin(), second.end(), intersection.begin());
        intersection.resize(a - intersection.begin());

        profile_diagonals.push_back(intersection.size());
    }
}


void GenMesh::stream_diagonals_out() {
    for (int j = 0; j < profile_diagonals.size(); j++)
        cout << j << ", "
             << radial_sorted_upper_bdr_indices[j] << ", " 
             << radial_sorted_lower_bdr_indices[j] << ", " 
             << profile_diagonals[j] << endl;
}


void GenMesh::stream_nodes_out() {
    ofstream nodes_file;
    nodes_file.open(output_dir + "/profile_nodes.dat", ios::trunc);
    for (auto& p: all_wafer_Point3D) {
        nodes_file << p.split(',') << endl;
    }
    nodes_file.close();
}


void GenMesh::stream_boundary_nodes_out() {

    ofstream sorted_boundary_nodes_file;
    sorted_boundary_nodes_file.open(output_dir + "/sorted_boundary_nodes_on_space.csv", ios::trunc);

    for (int l = number_of_layers - 1; l >= 0; l--) {

        for (auto& p: sorted_boundary_points) {
            sorted_boundary_nodes_file << Point3D(p.second, l * height_of_layer);
        }

    }
    sorted_boundary_nodes_file.close();
}


void GenMesh::construct_front_nodes_of_brick(vector<Point3D>& back_wall, int iteration) {

    /* back_wall represents a 4 x 3 matrix containing the coordinates
       the appended vector is a 4 x 3 matrix with the front coordinates

       First two are opposite to the midpoints
       Last  two are opposite to the two right boundary nodes: North East and South East.
    */

    Point3D middle_up   = (back_wall[3] + back_wall[0])/2;
    Point3D middle_dn   = (back_wall[2] + back_wall[1])/2;
    Point normal_dir_up = get_right_orthogonal(project2D(back_wall[3] - back_wall[0]));
    Point normal_dir_dn = get_right_orthogonal(project2D(back_wall[2] - back_wall[1]));
    Point front_point0  = project2D(middle_up) + normal_dir_up * delta_over_plane_xy[iteration];
    Point front_point1  = project2D(middle_dn) + normal_dir_dn * delta_over_plane_xy[iteration];
    Point radial_point0 = project2D(back_wall[3]) * (1 + delta_over_plane_xy[iteration]/project2D(back_wall[3]).norm());
    Point radial_point1 = project2D(back_wall[2]) * (1 + delta_over_plane_xy[iteration]/project2D(back_wall[2]).norm());
    vector<Point3D> front_nodes = {Point3D(front_point0,  upper_z[iteration]),
                                   Point3D(front_point1,  lower_z[iteration]),
                                   Point3D(radial_point0, upper_z[iteration]),
                                   Point3D(radial_point1, lower_z[iteration])};
    
    // here we add the actual coordinates of each new point.
    for (auto& node: front_nodes) {
        all_wafer_Point3D.push_back(node);
    }
}


void GenMesh::split_prism(int prisma_kind, vector<int> up, vector<int> dn) {
    /*
    
    This function decides the indices permutations to be added as new tetrahedra.

                    Diagonals on the back walls

                    up_0 ----- up_1 ----- up_2
    case 0           |     \     |   \      |      ==      (1, 1)           ----->         0
                     |      \    |    \     |
                    dn_0 ----- dn_1 ----- dn_2
    
                    up_0 ----- up_1 ----- up_2
    case 3           |      /    |   \      |      ==      (0, 1)           ----->         0
                     |     /     |    \     |      
                    dn_0 ----- dn_1 ----- dn_2



                    up_0 ----- up_1 ----- up_2
    case 1           |     /     |     /    |      ==      (0, 0)           ----->          2
                     |    /      |    /     |
                    dn_0 ----- dn_1 ----- dn_2

                    up_0 ----- up_1 ----- up_2
    case 2           |    \      |     /    |      ==      (1, 0)           ----->          2
                     |     \     |    /     |
                    dn_0 ----- dn_1 ----- dn_2

                    (false, false), (1, 0)                 (outer: false, true)
                    (false, true),  (3, 1)                 (outer: false, false)
                    (true, false),  (2, 0)                 (outer: false, true)
                    (true, true),   (0, 1)                 (outer: true, false)
    */

    vector<int> new_elem_0;
    vector<int> new_elem_1;
    vector<int> new_elem_2;

    switch (prisma_kind) {
        case 0:
            // The identity permutation.
            new_elem_0 = {up[0], up[1], up[2], dn[2]};
            new_elem_1 = {dn[0], dn[1], dn[2], up[0]};
            new_elem_2 = {up[0], up[1], dn[1], dn[2]};
            break;

        case 1:
            new_elem_0 = {up[0], up[1], up[2], dn[0]};
            new_elem_1 = {dn[0], dn[1], dn[2], up[2]};
            new_elem_2 = {dn[0], dn[1], up[1], up[2]};
            break;

        case 2:
            new_elem_0 = {up[0], up[1], up[2], dn[1]};
            new_elem_1 = {dn[0], dn[1], dn[2], up[2]};
            new_elem_2 = {up[0], up[2], dn[0], dn[1]};
            break;

        case 3:
            new_elem_0 = {up[0], up[1], up[2], dn[0]};
            new_elem_1 = {dn[0], dn[1], dn[2], up[1]};
            new_elem_2 = {up[1], up[2], dn[0], dn[2]};
            break;

        default:
            break;
    }

    // append the three new elements indices
    if (measure(new_elem_0) < 1e-12) {
        measure_flag = true;
        degenerate_elements.push_back(elements_by_vertices.size());
    }
    elements_by_vertices.emplace(elements_by_vertices.size(), new_elem_0);

    if (measure(new_elem_1) < 1e-12) {
        measure_flag = true;
        degenerate_elements.push_back(elements_by_vertices.size());
    }
    elements_by_vertices.emplace(elements_by_vertices.size(), new_elem_1);
    
    if (measure(new_elem_2) < 1e-12) {
        measure_flag = true;
        degenerate_elements.push_back(elements_by_vertices.size());
    }
    elements_by_vertices.emplace(elements_by_vertices.size(), new_elem_2);
}


void GenMesh::add_tetrahedra_within_prisms() {

    pair<bool, bool> current_inner_diagonals;
    pair<int, int> current_case;  // (left, right)
    pair<bool, bool> current_outer_diagonals;

    Prism left;
    Prism right;

    prisms_map.emplace(prisms_map.size(), prisms_map[0]);

    for (int j = 0; j < prisms_map.size() - 1; j++) {

        left         = prisms_map[j];
        right        = prisms_map[j + 1];
        
        current_inner_diagonals = make_pair(left.diagonal_orientation, 
                                            right.diagonal_orientation);

        current_case            = case_translator[current_inner_diagonals];
        current_outer_diagonals = outer_diagonals_translator[current_inner_diagonals];
        
        split_prism(current_case.first, 
                   {left.front_up, left.back[3], left.radial_up},
                   {left.front_dn, left.back[2], left.radial_dn});

        number_of_outer_walls++;
        next_profile_diagonals.push_back(current_outer_diagonals.first);
        new_indices_for_outer_facets.emplace(
            2 * j,
            vector<int>({left.front_up, left.front_dn, left.radial_dn, left.radial_up})
        );

        split_prism(current_case.second, 
                   {left.radial_up, left.back[3], right.front_up}, 
                   {left.radial_dn, left.back[2], right.front_dn});

        number_of_outer_walls++;
        next_profile_diagonals.push_back(current_outer_diagonals.second);
        new_indices_for_outer_facets.emplace(
            2 * j + 1,
            vector<int>({left.radial_up, left.radial_dn, right.front_dn, right.front_up})
        );
    }
}


void GenMesh::reset_profile_objects() {

    // erase
    radial_sorted_upper_bdr_indices = vector<int>(); 
    radial_sorted_lower_bdr_indices = vector<int>();

    // redefine
    indices_for_outer_facets = map<int, vector<int>>(new_indices_for_outer_facets);

    for (int b = 0; b < number_of_bricks_in_wall; b++) {
        radial_sorted_upper_bdr_indices.push_back(prisms_map[b].front_up);
        radial_sorted_upper_bdr_indices.push_back(prisms_map[b].radial_up);

        radial_sorted_lower_bdr_indices.push_back(prisms_map[b].front_dn);
        radial_sorted_lower_bdr_indices.push_back(prisms_map[b].radial_dn);
    }
    // Repeat first indices to have the last cuadrilateral.
    radial_sorted_upper_bdr_indices.push_back(radial_sorted_upper_bdr_indices[0]);
    radial_sorted_lower_bdr_indices.push_back(radial_sorted_lower_bdr_indices[0]);

    // erase
    new_indices_for_outer_facets    = map<int, vector<int>>();
    prisms_map                      = map<int, Prism>();
    sorted_boundary_points          = the_map();
    wall_of_Point3D                 = vector<Point3D>();
    physical_facets                 = map<int, vector<Point3D>>(); 
    
    number_of_bricks_in_wall        = number_of_outer_walls;
    profile_diagonals               = vector<bool>(next_profile_diagonals);

    next_profile_diagonals          = vector<bool>();
    number_of_outer_walls           = 0;

    for (int i = 0; i < number_of_bricks_in_wall; ++i) {
        physical_facets[i] = vector<Point3D>({
            all_wafer_Point3D[indices_for_outer_facets[i][0]],
            all_wafer_Point3D[indices_for_outer_facets[i][1]],
            all_wafer_Point3D[indices_for_outer_facets[i][2]],
            all_wafer_Point3D[indices_for_outer_facets[i][3]]
        });
    }

}

/*
void GenMesh::new_make_3D_points(vector<size_t>& in_vector) {
    double local_z;

    for (int l = number_of_layers - 1; l >= 0; l--) {
        local_z = l * height_of_layer;
        for (int j = 0; j < in_vector.size(); j++) {
            all_wafer_Point3D.push_back(Point3D(pointlist[in_vector[j]], local_z));
        }
    }
}

*/
void GenMesh::make_3D_points() {
    double local_z;

    for (int l = number_of_layers - 1; l >= 0; l--) {
        local_z = l * height_of_layer;
        for (int j = 0; j < pointlist.size(); j++) {
            all_wafer_Point3D.push_back(Point3D(pointlist[j], local_z));
        }
    }
}


void GenMesh::build_profile_mesh(int input_size) {

    for (int i = 0; i < input_size; i++)
    {
        upper_bdr_points_global_indices.push_back(i);
    }

    int added_index = all_wafer_Point3D.size();

    cout << "Making profile mesh starting with " << added_index << " boundary points." << endl;

    radial_sort_boundary_points();

    number_of_bricks_in_wall = sorted_boundary_points.size();
    bool this_diagonal       = false;
    double local_z;

    // Temporal objects to handle as current.
    vector<int> back_idxs(4);
    Prism current_prism;

    // Wall of 3D points for base case.
    // TODO: perhaps we can omit the construction of wall_of_Point3D
    // and do the base case physical_facets directly.
    for (int l = number_of_layers - 1; l >= 0; l--) {
        local_z = l * height_of_layer;
        for (auto& p: sorted_boundary_points) {
            wall_of_Point3D.push_back(Point3D(p.second, local_z));
        }
        wall_of_Point3D.push_back(Point3D(sorted_boundary_points[0], local_z));
    }

    for (int b = 0; b < number_of_bricks_in_wall; b++) {
        physical_facets[b] = vector<Point3D>({
            wall_of_Point3D[b],
            wall_of_Point3D[b + number_of_bricks_in_wall + 1],
            wall_of_Point3D[b + number_of_bricks_in_wall + 2],
            wall_of_Point3D[b + 1]
        });
    }
    
    // Iterate "control point levels"
    for (int c = 0; c < 5; c++) {

        if (debug) {
            cout << "CONTROL POINTS ITERATION NUMBER " << c + 1 << "\n";
            cout << "Current all_wafer_Point3D size " << added_index << endl;
            cout << "Current element number " << elements_by_vertices.size() << endl;
            cout << "number_of_bricks_in_wall " << number_of_bricks_in_wall << "\n";
            cout << "physical_facets " << physical_facets.size() << "\n";
            cout << "profile_diagonals " << profile_diagonals.size() << "\n";
            cout << "Sizes of radial--sorted indices (" 
                 << radial_sorted_upper_bdr_indices.size() << ", "
                 << radial_sorted_lower_bdr_indices.size() << ") " << endl;
        }

        for (int b = 0; b < number_of_bricks_in_wall; b++) {
            // Append actual R3 points to nodes.csv
            construct_front_nodes_of_brick(physical_facets[b], c);
            this_diagonal = profile_diagonals[b];
            back_idxs = {radial_sorted_upper_bdr_indices[b],
                         radial_sorted_lower_bdr_indices[b],
                         radial_sorted_lower_bdr_indices[b + 1],
                         radial_sorted_upper_bdr_indices[b + 1]};

            current_prism.diagonal_orientation = this_diagonal;
            current_prism.front_up             = added_index;
            current_prism.front_dn             = added_index + 1;
            current_prism.back                 = back_idxs;

            // Indices of the points opposite to the boundary points
            current_prism.radial_up            = added_index + 2;
            current_prism.radial_dn            = added_index + 3;

            prisms_map.emplace(b, current_prism);

            // mesh the first array of prisms
            if (this_diagonal) {
                split_prism(2, {back_idxs[0], back_idxs[3], current_prism.front_up},
                               {back_idxs[1], back_idxs[2], current_prism.front_dn});
            } else {
                split_prism(1, {current_prism.front_up, back_idxs[0], back_idxs[3]},
                               {current_prism.front_dn, back_idxs[1], back_idxs[2]});
            }
            added_index += 4;

        }

        add_tetrahedra_within_prisms();
        reset_profile_objects();

    }

    cout << "Done.\n\n";

}

// TODO: after Type 2
// double max()
// double min()
// report the max of the whole mesh
// report the min for the whole ...

double GenMesh::measure(vector<int>& tetra) {
    Point3D i_versor = all_wafer_Point3D[tetra[1]] - all_wafer_Point3D[tetra[0]];
    Point3D j_versor = all_wafer_Point3D[tetra[2]] - all_wafer_Point3D[tetra[0]];
    Point3D k_versor = all_wafer_Point3D[tetra[3]] - all_wafer_Point3D[tetra[0]];

    return fabs(i_versor * j_versor.product(k_versor))/6;
}
