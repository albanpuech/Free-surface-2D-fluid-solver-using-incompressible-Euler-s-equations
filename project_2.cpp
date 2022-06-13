#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <string>
#include <iostream>
#include <initializer_list>
#include <string>
#include <algorithm>
#include <tuple>
#include <list>
#include <random>
#include <chrono>
#include "lbfgs.h"
#include <sstream>

int sgn(double value)
{
    if (value > 0)
        return 1;
    if (value < 0)
        return -1;
    return 0;
}

// vector class from previous assignment
class Vector
{
public:
    explicit Vector(double x = 0, double y = 0)
    {
        data[0] = x;
        data[1] = y;
    }
    double norm2() const
    {
        return data[0] * data[0] + data[1] * data[1];
    }
    double norm() const
    {
        return sqrt(norm2());
    }
    void normalize()
    {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
    }
    double operator[](int i) const { return data[i]; };
    double &operator[](int i) { return data[i]; };
    double data[2];
};

Vector operator+(const Vector &a, const Vector &b)
{
    return Vector(a[0] + b[0], a[1] + b[1]);
}
Vector operator-(const Vector &a, const Vector &b)
{
    return Vector(a[0] - b[0], a[1] - b[1]);
}
Vector operator-(const Vector &a)
{
    return Vector(-a[0], -a[1]);
}
Vector operator*(const double a, const Vector &b)
{
    return Vector(a * b[0], a * b[1]);
}
Vector operator*(const Vector &a, const double b)
{
    return Vector(a[0] * b, a[1] * b);
}
Vector operator/(const Vector &a, const double b)
{
    return Vector(a[0] / b, a[1] / b);
}

double dot(const Vector &a, const Vector &b)
{

    return a[0] * b[0] + a[1] * b[1];
}

Vector operator*(const Vector &a, const Vector &b)
{
    return Vector(a[0] * b[0], a[1] * b[1]);
}

// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name
class Polygon
{
public:
    double area()
    {
        int N = vertices.size();
        if (N <= 2)
            return 0;
        double area_val = 0.;
        for (int i = 0; i < N; ++i)
        {
            area_val += ((vertices[i][0] * vertices[(i + 1) % N][1]) - (vertices[i][1] * vertices[(i + 1) % N][0]));
        }
        return std::abs(area_val / 2);
    }

    double integral_squared_dist(const Vector &P)
    {
        double dist = 0;
        if (vertices.size() <= 2)
        {
            return dist;
        }
        for (int i = 1; i < vertices.size() - 1; ++i)
        {
            double area = 0.5 * abs((vertices[i][0] - vertices[0][0]) * (vertices[i + 1][1] - vertices[0][1]) - (vertices[i][1] - vertices[0][1]) * (vertices[i + 1][0] - vertices[0][0]));
            Vector C[3] = {vertices[0], vertices[i], vertices[i + 1]};
            for (int k = 0; k < 3; ++k)
            {
                for (int l = k; l < 3; ++l)
                {
                    dist += area / 6 * dot(C[k] - P, C[l] - P);
                }
            }
        }
        return std::abs(dist);
    }

    Vector centroid()
    {
        int N = vertices.size();
        Vector Center(0, 0);
        double area_val = area();
        if (area_val == 0)
        {
            if (N != 0)
            {
                return vertices[0];
            }
            else
            {
                return Vector(0, 0);
            }
        }
        for (int i = 0; i < N; i++)
        {
            Center[0] += (vertices[i][0] + vertices[(i + 1) % N][0]) * (vertices[i][0] * vertices[(i + 1) % N][1] - vertices[(i + 1) % N][0] * vertices[i][1]);
            Center[1] += (vertices[i][1] + vertices[(i + 1) % N][1]) * (vertices[i][0] * vertices[(i + 1) % N][1] - vertices[(i + 1) % N][0] * vertices[i][1]);
        }
        return -Center / (6. * area_val);
    }
    std::vector<Vector> vertices;
};

// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none")
{
    FILE *f = fopen(filename.c_str(), "w+");
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (int i = 0; i < polygons.size(); i++)
    {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \"");
        for (int j = 0; j < polygons[i].vertices.size(); j++)
        {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}

Polygon clip_cell(const Polygon &cell, int i, int j, const Vector *points, const double *weights, int N)
{
    Polygon res_poly;
    for (int k = 0; k < cell.vertices.size(); ++k)
    {
        const Vector &A = cell.vertices[((k - 1) >= 0) ? (k - 1) : (cell.vertices.size() - 1)];
        const Vector &B = cell.vertices[k];

        Vector M = (points[i] + points[j]) / 2 + (weights[i] - weights[j]) / (2 * (points[i] - points[j]).norm2()) * (points[j] - points[i]);

        double t = dot(M - A, points[j] - points[i]) / dot(B - A, points[j] - points[i]);
        Vector P = A + t * (B - A);

        if ((B - points[i]).norm2() - weights[i] <= (B - points[j]).norm2() - weights[j])
        {
            if ((A - points[i]).norm2() - weights[i] > (A - points[j]).norm2() - weights[j])
            {
                res_poly.vertices.push_back(P);
            }
            res_poly.vertices.push_back(B);
        }
        else
        {
            if ((A - points[i]).norm2() - weights[i] <= (A - points[j]).norm2() - weights[j])
            {
                res_poly.vertices.push_back(P);
            }
        }
    }
    return res_poly;
}

Polygon clip_polygon_by_edge(const Polygon &cell, const Vector &poly1, const Vector &poly2)
{
    Polygon res_poly;
    Vector u = poly1;
    Vector n(poly2[1] - poly1[1], -poly2[0] + poly1[0]);
    int N = cell.vertices.size();

    for (int k = 0; k < N; ++k)
    {
        const Vector &A = cell.vertices[((k - 1) >= 0) ? (k - 1) : (N - 1)];
        const Vector &B = cell.vertices[k];

        double t = dot(u - A, n) / dot(B - A, n);
        Vector P = A + t * (B - A);

        if (dot(u - B, n) <= 0)
        {
            if (dot(u - A, n) > 0)
            {
                res_poly.vertices.push_back(P);
            }
            res_poly.vertices.push_back(B);
        }
        else
        {
            if (dot(u - A, n) <= 0)
            {
                res_poly.vertices.push_back(P);
            }
        }
    }
    return res_poly;
}

Polygon clip_poly(const Polygon &poly, const Polygon &clip)
{
    Polygon result = poly;
    for (int i = 0; i < clip.vertices.size(); i++)
    {
        result = clip_polygon_by_edge(result, clip.vertices[i], clip.vertices[(i + 1) % clip.vertices.size()]);
    }
    return result;
}

Polygon compute_cell(int i, const Vector *points, const double *weights, int N)
{
    Polygon cell;
    cell.vertices.push_back(Vector(0, 0));
    cell.vertices.push_back(Vector(0, 1));
    cell.vertices.push_back(Vector(1, 1));
    cell.vertices.push_back(Vector(1, 0));

    for (int j = 0; j < N; ++j)
    {
        cell = clip_cell(cell, i, j, points, weights, N);
    }
    return cell;
}

std::vector<Polygon> get_diagram(const Vector *points, const double *w, int N)
{
    std::vector<Polygon> diagram(N);
    for (int i = 0; i < N; ++i)
    {
        diagram[i] = compute_cell(i, &points[0], &w[0], N);
    }
    return diagram;
}

std::vector<Polygon> get_diagram_fluid(const Vector *points, const double *w, int N)
{
    std::vector<Polygon> diagram(N);
#pragma omp parallel for
    for (int i = 0; i < N - 1; ++i)
    {

        Polygon disc_poly;
        disc_poly.vertices.resize(100);
        double radius = sqrt(w[i] - w[N - 1]);
        for (int j = 0; j < 100; j++)
        {
            disc_poly.vertices[j][0] = cos(j / 100. * 2 * M_PI) * radius + points[i][0];
            disc_poly.vertices[j][1] = -sin(j / 100. * 2 * M_PI) * radius + points[i][1];
        }
        diagram[i] = compute_cell(i, &points[0], &w[0], N - 1);
        diagram[i] = clip_poly(diagram[i], disc_poly);
    }
    return diagram;
}

lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step)
{
    int i;
    lbfgsfloatval_t fx = 0.0;

    Vector *points = static_cast<Vector *>(instance);
    std::vector<Polygon> diagram = get_diagram(&points[0], &x[0], n);

    double lambda = 1.0 / n;
    for (int i = 0; i < n; i++)
    {
        double area = diagram[i].area();
        g[i] = -(lambda - area);
        fx += -(diagram[i].integral_squared_dist(points[i]) - x[i] * area + lambda * x[i]);
    }

    return fx;
}

lbfgsfloatval_t evaluate_fluid(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step)
{
    int i;
    lbfgsfloatval_t fx = 0.0;

    Vector *points = static_cast<Vector *>(instance);
    std::vector<Polygon> diagram = get_diagram_fluid(&points[0], &x[0], n);
    double fraction_fluid = 0.4;
    double fraction_air = 1 - fraction_fluid;
    double lambda = fraction_fluid / (n - 1);
    double sum_fluid_area = 0;
    for (int i = 0; i < n - 1; i++)
    {
        double area = diagram[i].area();
        sum_fluid_area += area;
        g[i] = -(lambda - area);
        fx += -(diagram[i].integral_squared_dist(points[i]) - x[i] * area + lambda * x[i]);
    }
    // last component
    double air_area = 1 - sum_fluid_area;
    g[n - 1] = -(fraction_air - air_area);
    fx += -(-x[n - 1] * air_area + fraction_air * x[n - 1]);

    return fx;
}

void fluid_time_step(std::vector<Polygon> &diagram, std::vector<Vector> &particles, std::vector<Vector> &vel, std::vector<double> &w)
{

    double fx;

    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;

    int N = particles.size() + 1;

    int ret = lbfgs(N, &w[0], &fx, evaluate_fluid, nullptr, &particles[0], &param);
    diagram = get_diagram_fluid(&particles[0], &w[0], N);
    // moving points is now done in main function so that it occurs after we optimize the weights and compute the new power diagram
}
void save_frame(const std::vector<Polygon> &cells, std::string filename, int frameid = 0)
{
    int W = 1000, H = 1000;
    std::vector<unsigned char> image(W * H * 3, 255);
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < cells.size(); i++)
    {

        double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;   
        for (int j = 0; j < cells[i].vertices.size(); j++)
        {
            bminx = std::min(bminx, cells[i].vertices[j][0]);
            bminy = std::min(bminy, cells[i].vertices[j][1]);
            bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
            bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
        }
        bminx = std::min(W - 1., std::max(0., W * bminx));
        bminy = std::min(H - 1., std::max(0., H * bminy));
        bmaxx = std::max(W - 1., std::max(0., W * bmaxx));
        bmaxy = std::max(H - 1., std::max(0., H * bmaxy));

        for (int y = bminy; y < bmaxy; y++)
        {
            for (int x = bminx; x < bmaxx; x++)
            {
                int prevSign = 0;
                bool isInside = true;
                double mindistEdge = 1E9;
                for (int j = 0; j < cells[i].vertices.size(); j++)
                {
                    double x0 = cells[i].vertices[j][0] * W;
                    double y0 = cells[i].vertices[j][1] * H;
                    double x1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][0] * W;
                    double y1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][1] * H;
                    double det = (x - x0) * (y1 - y0) - (y - y0) * (x1 - x0);
                    int sign = sgn(det);
                    if (prevSign == 0)
                        prevSign = sign;
                    else if (sign == 0)
                        sign = prevSign;
                    else if (sign != prevSign)
                    {
                        isInside = false;
                        break;
                    }
                    prevSign = sign;
                    double edgeLen = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
                    double distEdge = std::abs(det) / edgeLen;
                    double dotp = (x - x0) * (x1 - x0) + (y - y0) * (y1 - y0);
                    if (dotp < 0 || dotp > edgeLen * edgeLen)
                        distEdge = 1E9;
                    mindistEdge = std::min(mindistEdge, distEdge);
                }
                if (isInside)
                {
                    if (i < 200)
                    { // the N first particles may represent fluid, displayed in blue
                        image[((H - y - 1) * W + x) * 3] = 0;
                        image[((H - y - 1) * W + x) * 3 + 1] = 0;
                        image[((H - y - 1) * W + x) * 3 + 2] = 255;
                    }
                    if (mindistEdge <= 2)
                    {
                        image[((H - y - 1) * W + x) * 3] = 0;
                        image[((H - y - 1) * W + x) * 3 + 1] = 0;
                        image[((H - y - 1) * W + x) * 3 + 2] = 0;
                    }
                }
            }
        }
    }
    std::ostringstream os;
    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    int N = 200;
    std::vector<Vector> particles(N);
    std::vector<Vector> vel(N, Vector(0, 0));
    std::vector<Polygon> diagram;

    double m = 200;
    double eps2 = 0.004 * 0.004;
    double dt = 0.004;

    std::vector<double> w(N);
    for (int i = 0; i < N; ++i)
    {
        particles[i][0] = rand() / (double)RAND_MAX;
        particles[i][1] = rand() / (double)RAND_MAX;
        w[i] = 0.005;
    }

    for (int t = 0; t < 500; t++)
    {
        // optimize weights
        fluid_time_step(diagram, particles, vel, w);
        // compute diagram
        std::vector<Polygon> cells = get_diagram_fluid(&particles[0], &w[0], N + 1);
        // save frame
        save_frame(cells, "frame_", t);
        // move points
        for (int i = 0; i < particles.size(); ++i)
        {
            Vector spring_force = 1 / eps2 * (diagram[i].centroid() - particles[i]);
            Vector forces = spring_force + m * Vector(0, -9.81);
            vel[i] = vel[i] + dt / m * forces;
            particles[i] = particles[i] + dt * vel[i];
            if (particles[i][0] < 0)
            {
                particles[i][0] = -particles[i][0];
                vel[i][0] = 0;
            }
            if (particles[i][1] < 0)
            {
                particles[i][1] = -particles[i][1];
                vel[i][1] = 0;
            }
            if (particles[i][0] > 1)
            {
                particles[i][0] = 2 - particles[i][0];
                vel[i][0] = 0;
            }
            if (particles[i][1] > 1)
            {
                particles[i][1] = 2 - particles[i][1];
                vel[i][1] = 0;
            }
        }
    }
        auto end = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Rendering Time: " << time.count() << "ms" << std::endl;
    }

    // 	std::vector<Polygon> diagram = get_diagram(&points[0], &w[0], N);
    // 	save_svg(diagram, "diagram_before.svg");
    // 	// #pragma omp parallel for schedule(dynamic, 1)

    // 	lbfgs_parameter_t param;
    // 	lbfgs_parameter_init(&param);
    // 	param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
    // 	double fx;
    // 	int ret = lbfgs(N, &w[0], &fx, evaluate, progress, &points[0], &param);
    // 	std::cout << ret << std::endl;
    // 	diagram = get_diagram(&points[0], &w[0], N);

    // 	save_svg(diagram, "diagram.svg");

    //     auto end = std::chrono::high_resolution_clock::now();
    //     auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    //     std::cout << "Rendering Time: " << time.count() << "ms" << std::endl;

    // 	return 0;
    // }
