#include "GenMesh.h"


GenMesh::GenMesh(
    int n_layers,
    double h_layer,
    bool debug_flag,
    string prefix,
    map<string, valarray<double>>& profile_params
):

    number_of_layers(n_layers),
    height_of_layer(h_layer),
    debug(debug_flag),
    output_dir(prefix)

{

    // TODO: remove these three, as we can call the map<> profile_params directly
    delta_over_plane_xy = profile_params["Width"];
    upper_z             = profile_params["Ceiling"];
    lower_z             = profile_params["Floor"];

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


void find_global_coordinates_for_boundary()
{
    
}


Point GenMesh::get_right_orthogonal(Point vector2D)
{
    // TODO: test once more
    return Point(vector2D.m_y, -vector2D.m_x)/vector2D.norm();
}


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

// TODO: to study the condition of the matrix: 
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
