#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cstring>
#include <string>
#include <windows.h>
using namespace std;

# define PI           3.14159265358979323846

struct Vec3
{
    float x, y, z;

    Vec3(float x, float y, float z): x(x), y(y), z(z)
    {}
    
    Vec3 operator + (const Vec3& other)
    {
        return Vec3(x + other.x, y + other.y, z + other.z);
    }

    Vec3 operator - (const Vec3& other)
    {
        return Vec3(x - other.x, y - other.y, z - other.z);
    }

    float dot(const Vec3& other)
    {
        return x * other.x + y * other.y + z * other.z;
    }

    Vec3 normalize()
    {
        float magnitude = sqrt(x*x + y*y + z*z);
        if(magnitude == 0.0f) return Vec3(0,0,0);
        return Vec3(x / magnitude, y / magnitude, z / magnitude);
    }

};

struct Mat3
{
    float m[3][3];


    Vec3 operator * (const Vec3& v)
    {
        return 
        {
            m[0][0]*v.x + m[0][1]*v.y + m[0][2]*v.z,
            m[1][0]*v.x + m[1][1]*v.y + m[1][2]*v.z,
            m[2][0]*v.x + m[2][1]*v.y + m[2][2]*v.z
        };
    }

    Mat3 operator * (const Mat3& other)
    {
        Mat3 result{};
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                result.m[i][j] = 0;
                for (int k = 0; k < 3; ++k)
                {
                    result.m[i][j] += m[i][k] * other.m[k][j];
                }
            }
        }
        return result;
    }
};

// constants
const int WIDTH = 80, HEIGHT = 40;
const float R1 = 0.5; // small radius
const float R2 = 1.0; // large radius
const float K2 = 10.0; // distance from camera
const float K1 = (WIDTH * K2 * 3 / (8 * (R1 + R2))) * 0.6; // scaling factor
const string brightness = ".,-~:;=!*#$@";
// light direction vector
const Vec3 lightDirection  = Vec3(1,1,-2).normalize();

Mat3 rotationMatrixX(float angle)
{
    float c = cos(angle);
    float s = sin(angle);
    Mat3 mat;
    mat.m[0][0] = 1; mat.m[0][1] = 0; mat.m[0][2] = 0;
    mat.m[1][0] = 0; mat.m[1][1] = c; mat.m[1][2] = -s;
    mat.m[2][0] = 0; mat.m[2][1] = s; mat.m[2][2] = c;
    return mat;
}

Mat3 rotationMatrixY(float angle)
{
    float c = cos(angle);
    float s = sin(angle);
    Mat3 mat;
    mat.m[0][0] = c; mat.m[0][1] = 0; mat.m[0][2] = s;
    mat.m[1][0] = 0; mat.m[1][1] = 1; mat.m[1][2] = 0;
    mat.m[2][0] = -s; mat.m[2][1] = 0; mat.m[2][2] = c;
    return mat;
}

Mat3 rotationMatrixZ(float angle)
{
    float c = cos(angle);
    float s = sin(angle);
    Mat3 mat;
    mat.m[0][0] = c; mat.m[0][1] = -s; mat.m[0][2] = 0;
    mat.m[1][0] = s; mat.m[1][1] = c; mat.m[1][2] = 0;
    mat.m[2][0] = 0; mat.m[2][1] = 0; mat.m[2][2] = 1;
    return mat;
}

void render(float A, float B)
{
    vector<char> output(WIDTH * HEIGHT , ' ');
    vector<float> zBuffer(WIDTH * HEIGHT, 0);

    Mat3 rotationX = rotationMatrixX(A);
    Mat3 rotationZ = rotationMatrixZ(B);
    Mat3 rotationY = rotationMatrixY(B);
    Mat3 combinedRotation = rotationX * rotationY;

    for (float theta = 0; theta < 2 * PI; theta += 0.1f) // 0.07 originally
    {
        float cosTheta = cos(theta), sinTheta = sin(theta);
        for (float delta = 0; delta < 2 * PI; delta += 0.05f) // 0.02 originally
        {
            float cosDelta = cos(delta), sinDelta = sin(delta);
            // find point on donut
            Vec3 center(R2 * cosDelta, R2 * sinDelta, 0); // point on the circle  
            Vec3 outward(R1 * cosTheta * cosDelta, R1 * cosTheta * sinDelta, R1 * sinTheta);

            Vec3 surfacePoint = center + outward;
            Vec3 normal(cosTheta * cosDelta, cosTheta * sinDelta, sinTheta);

            Vec3 rotatedSurfacePoint = combinedRotation * surfacePoint;
            Vec3 rotatedNormal = combinedRotation * normal;
            
            float luminance = rotatedNormal.normalize().dot(lightDirection);

            if (luminance > 0)
            {
                // project 3D point onto x-y plane with perpective
                // where K2 is distanc from camera // over 1 ro simulate perpective
                float projZ = 1.0 / (K2 + rotatedSurfacePoint.z);

                // where K1 is used to scale point to fit in screen
                int projX = static_cast<int>(WIDTH / 2 + K1 * projZ * rotatedSurfacePoint.x);
                int projY = static_cast<int>(HEIGHT / 2 - K1 * projZ * rotatedSurfacePoint.y);

                int bufferIndex = projX + projY * WIDTH;

                if (bufferIndex >= 0 && bufferIndex < WIDTH * HEIGHT && projZ > zBuffer[bufferIndex])
                {
                    zBuffer[bufferIndex] = projZ;

                    // map luminance to ascii
                    int brightnessIndex = (int)(luminance * (brightness.length() - 1));
                    brightnessIndex = max(0, min((int)brightness.length() - 1, brightnessIndex));
                    output[bufferIndex] = brightness[brightnessIndex];
                }

            }
        }
    }

    cout << "\x1b[H"; // reset cursor

    for (int j = 0; j < HEIGHT; ++j)
    {
        for (int i = 0; i < WIDTH; ++i)
        {
            cout << output[i + j * WIDTH];
        }
        cout << endl;
    }

}

int main() 
{
    // chatgpt code to reset cursor
    HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
    if (hOut == INVALID_HANDLE_VALUE) {
        return GetLastError(); // Handle error
    }
    DWORD dwMode = 0;
    if (!GetConsoleMode(hOut, &dwMode)) {
        return GetLastError(); // Handle error
    }
    dwMode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
    if (!SetConsoleMode(hOut, dwMode)) {
        // This might fail on older Windows versions, proceed anyway
        // or fallback to another method if necessary.
        // return GetLastError(); // Optional: Handle error strictly
    }
    // end of reset cursor code

    float A = 0, B = 0; // Rotation angles

    cout << "Starting Donut..." << endl;
    Sleep(1000); // pause 1 sec

    // animation loop
    while(true)
    {
        render(A, B);
        A += 0.04f; // rotation angle A
        B += 0.02f; // rotation angle B

        // add delay
        //Sleep(1); // Adjust for desired speed (microseconds)
        
    }

    return 0;
}
