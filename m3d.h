//
// Created by genkinger on 9/27/17.
// VERSION 0.1
//

#ifndef M3D_M3D_H
#define M3D_M3D_H

#include <math.h>
#include <stdio.h>


typedef struct{
    union {int x,u;};
    union {int y,v;};
} vec2i_t;

typedef struct{
    union {int x,r;};
    union {int y,g;};
    union {int z,b;};
} vec3i_t;

typedef struct{
    union {int x,r;};
    union {int y,g;};
    union {int z,b;};
    union {int w,a;};
} vec4i_t;


typedef struct{
    union {float x, u;};
    union {float y, v;};
} vec2f_t;

typedef struct{
    union {float x,r;};
    union {float y,g;};
    union {float z,b;};
} vec3f_t;

typedef struct{
    union {float x,r;};
    union {float y,g;};
    union {float z,b;};
    union {float w,a;};
} vec4f_t;

typedef union{
        float elements[4][4];
        struct{
            float m00,m01,m02,m03;
            float m10,m11,m12,m13;
            float m20,m21,m22,m23;
            float m30,m31,m32,m33;
        };
} mat4_t;

typedef union{
    float elements[3][3];
    struct {
        float m00,m01,m02;
        float m10,m11,m12;
        float m20,m21,m22;
    };
} mat3_t;

#ifdef M3D_ENABLE_INTEGER_FUNCTIONS
    vec2i_t vec2i(int x, int y);
    vec3i_t vec3i(int x, int y, int z);
    vec4i_t vec4i(int x, int y, int z, int w);
    //TODO: add all other functions here !!
#endif

mat4_t mat4_identity(float diagonal);
mat3_t mat3_identity(float diagonal);

mat4_t mat4(float m00, float m01, float m02, float m03,
            float m10, float m11, float m12, float m13,
            float m20, float m21, float m22, float m23,
            float m30, float m31, float m32, float m33);
mat3_t mat3(float m00, float m01, float m02,
            float m10, float m11, float m12,
            float m20, float m21, float m22);

vec2f_t vec2f(float x, float y);
vec3f_t vec3f(float x, float y, float z);
vec4f_t vec4f(float x, float y, float z, float w);




/**** OPERATIONS  vec4f ****/
vec4f_t vec4f_multiply_mat4(vec4f_t vec, mat4_t mat);
float   vec4f_dot_vec4f(vec4f_t first, vec4f_t second);
float   vec4f_length(vec4f_t vec);
vec4f_t vec4f_normalize_copy(vec4f_t vec);
void    vec4f_normalize(vec4f_t* vec);
vec4f_t vec4f_add_vec4f(vec4f_t first, vec4f_t second);
vec4f_t vec4f_subtract_vec4f(vec4f_t first, vec4f_t second);
vec4f_t vec4f_multiply_vec4f(vec4f_t first, vec4f_t second);
vec4f_t vec4f_divide_vec4f(vec4f_t first, vec4f_t second);
vec4f_t vec4f_add_float(vec4f_t vec, float value);
vec4f_t vec4f_subtract_float(vec4f_t vec, float value);
vec4f_t vec4f_multiply_float(vec4f_t vec, float value);
vec4f_t vec4f_divide_float(vec4f_t vec, float value);


/**** OPERATIONS  vec3f ****/
vec3f_t vec3f_multiply_mat3(vec3f_t vec, mat3_t mat);
float   vec3f_dot_vec3f(vec3f_t first, vec3f_t second);
vec3f_t   vec3f_cross_vec3f(vec3f_t first, vec3f_t second);
float   vec3f_length(vec3f_t vec);
vec3f_t vec3f_normalize_copy(vec3f_t vec);
void    vec3f_normalize(vec3f_t* vec);
vec3f_t vec3f_add_vec3f(vec3f_t first, vec3f_t second);
vec3f_t vec3f_subtract_vec3f(vec3f_t first, vec3f_t second);
vec3f_t vec3f_multiply_vec3f(vec3f_t first, vec3f_t second);
vec3f_t vec3f_divide_vec3f(vec3f_t first, vec3f_t second);
vec3f_t vec3f_add_float(vec3f_t vec, float value);
vec3f_t vec3f_subtract_float(vec3f_t vec, float value);
vec3f_t vec3f_multiply_float(vec3f_t vec, float value);
vec3f_t vec3f_divide_float(vec3f_t vec, float value);


/**** OPERATIONS  vec2f ****/
float   vec2f_dot_vec2f(vec2f_t first, vec2f_t second);
float   vec2f_length(vec2f_t vec);
vec2f_t vec2f_normalize_copy(vec2f_t vec);
void    vec2f_normalize(vec2f_t* vec);
vec2f_t vec2f_add_vec2f(vec2f_t first, vec2f_t second);
vec2f_t vec2f_subtract_vec2f(vec2f_t first, vec2f_t second);
vec2f_t vec2f_multiply_vec2f(vec2f_t first, vec2f_t second);
vec2f_t vec2f_divide_vec2f(vec2f_t first, vec2f_t second);
vec2f_t vec2f_add_float(vec2f_t vec, float value);
vec2f_t vec2f_subtract_float(vec2f_t vec, float value);
vec2f_t vec2f_multiply_float(vec2f_t vec, float value);
vec2f_t vec2f_divide_float(vec2f_t vec, float value);


/**** OPERATIONS mat4 ****/
mat4_t mat4_multiply_mat4(mat4_t first, mat4_t second);
mat4_t mat4_invert_affine(mat4_t matrix);
mat4_t mat4_transpose(mat4_t matrix);
mat4_t mat4_perspective_deg(float fovy, float near, float far, float aspect);
mat4_t mat4_perspective_rad(float fovy, float near, float far, float aspect);
mat4_t mat4_look_at(vec3f_t from, vec3f_t to, vec3f_t up);
mat4_t mat4_orthogonal(float left, float right, float top, float bottom, float back, float front);
mat4_t mat4_rotation_deg(float angle,vec3f_t axis);
mat4_t mat4_rotation_rad(float angle,vec3f_t axis);
mat4_t mat4_translation(vec3f_t translation);
mat4_t mat4_scaling(vec3f_t scale);

int    m3d_in_range(float value, float min, float max);
float  m3d_map(float value, float in_min, float in_max, float out_min, float out_max);
float  m3d_constrain(float value, float min, float max);
float  m3d_interpolate_linear(float start, float end, float factor);

/**** Debug ****/
void m3d_print_vec2f_raw(vec2f_t vec, const char* name);
void m3d_print_vec3f_raw(vec3f_t vec, const char* name);
void m3d_print_vec4f_raw(vec4f_t vec, const char* name);
void m3d_print_mat4_raw(mat4_t mat, const char* name);
void m3d_print_mat3_raw(mat3_t mat, const char* name);

#define m3d_print_vec2f(arg) m3d_print_vec2f_raw(arg,#arg)
#define m3d_print_vec3f(arg) m3d_print_vec3f_raw(arg,#arg)
#define m3d_print_vec4f(arg) m3d_print_vec4f_raw(arg,#arg)
#define m3d_print_mat4(arg) m3d_print_mat4_raw(arg,#arg)
#define m3d_print_mat3(arg) m3d_print_mat3_raw(arg,#arg)



/**** handy defines *****/
#define MAT4() mat4_identity(1.0f)
#define MAT3() mat3_identity(1.0f)
#define VEC3F() vec3f(0,0,0)
#define VEC4F() vec4f(0,0,0,0)
#define VEC2F() vec2f(0,0)
#define M3D_PI 3.1415926535897932384626433


#define DEG_TO_RAD(x)  x * M3D_PI/180.0
#define RAD_TO_DEG(x)  x * 180.0/M3D_PI

#define MAX(a,b) ((a > b) ? a : b)
#define MIN(a,b) ((a < b) ? a : b)


#endif //M3D_M3D_H
