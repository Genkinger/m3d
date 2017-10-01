//
// Created by genkinger on 9/27/17.
// VERSION: 0.1
//


#include "m3d.h"

mat4_t mat4_identity(float diagonal){
    mat4_t m;
    for (int i = 0; i < 4; ++i) {
        m.elements[i][0] = 0;
        m.elements[i][1] = 0;
        m.elements[i][2] = 0;
        m.elements[i][3] = 0;
    }
    m.elements[0][0] = diagonal;
    m.elements[1][1] = diagonal;
    m.elements[2][2] = diagonal;
    m.elements[3][3] = diagonal;
    return m;
}
mat3_t mat3_identity(float diagonal){
    mat3_t m;
    for (int i = 0; i < 3; ++i) {
        m.elements[i][0] = 0;
        m.elements[i][1] = 0;
        m.elements[i][2] = 0;
    }
    m.elements[0][0] = diagonal;
    m.elements[1][1] = diagonal;
    m.elements[2][2] = diagonal;
    return m;
}
mat4_t mat4(float m00, float m01, float m02, float m03,
            float m10, float m11, float m12, float m13,
            float m20, float m21, float m22, float m23,
            float m30, float m31, float m32, float m33)
{
    mat4_t mat;
    mat.elements[0][0] = m00;
    mat.elements[0][1] = m01;
    mat.elements[0][2] = m02;
    mat.elements[0][3] = m03;

    mat.elements[1][0] = m10;
    mat.elements[1][1] = m11;
    mat.elements[1][2] = m12;
    mat.elements[1][3] = m13;

    mat.elements[2][0] = m20;
    mat.elements[2][1] = m21;
    mat.elements[2][2] = m22;
    mat.elements[2][3] = m23;

    mat.elements[3][0] = m30;
    mat.elements[3][1] = m31;
    mat.elements[3][2] = m32;
    mat.elements[3][3] = m33;

    return mat;
}
mat3_t mat3(float m00, float m01, float m02,
            float m10, float m11, float m12,
            float m20, float m21, float m22)
{
    mat3_t mat;
    mat.elements[0][0] = m00;
    mat.elements[0][1] = m01;
    mat.elements[0][2] = m02;


    mat.elements[1][0] = m10;
    mat.elements[1][1] = m11;
    mat.elements[1][2] = m12;


    mat.elements[2][0] = m20;
    mat.elements[2][1] = m21;
    mat.elements[2][2] = m22;
    return mat;
}


vec2f_t vec2f(float x, float y){
    return (vec2f_t){x,y};
}
vec3f_t vec3f(float x, float y, float z){
    return (vec3f_t){x,y,z};
}
vec4f_t vec4f(float x, float y, float z, float w){
    return (vec4f_t){x,y,z,w};
}



/**** OPERATIONS  vec4f ****/
vec4f_t vec4f_multiply_mat4(vec4f_t vec, mat4_t mat){
    vec4f_t a = vec4f_multiply_float(vec4f(mat.m00,mat.m10,mat.m20,mat.m30),vec.x);
    vec4f_t b = vec4f_multiply_float(vec4f(mat.m01,mat.m11,mat.m21,mat.m31),vec.y);
    vec4f_t c = vec4f_multiply_float(vec4f(mat.m02,mat.m12,mat.m22,mat.m32),vec.z);
    vec4f_t d = vec4f_multiply_float(vec4f(mat.m03,mat.m13,mat.m23,mat.m33),vec.w);

    return vec4f_add_vec4f(a,vec4f_add_vec4f(b,vec4f_add_vec4f(c,d)));

}
float   vec4f_dot_vec4f(vec4f_t first, vec4f_t second){
    return (first.x * second.x + first.y * second.y + first.z * second.z + first.w * second.w);
}
float   vec4f_length(vec4f_t vec){
    return sqrtf(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z + vec.w*vec.w);
}
vec4f_t vec4f_normalize_copy(vec4f_t vec){
    float l = vec4f_length(vec);
    if(l > 0)
        return vec4f(vec.x/l,vec.y/l,vec.z/l,vec.w/l);
    else
        return vec4f(0,0,0,0);
}
void    vec4f_normalize(vec4f_t* vec){}
vec4f_t vec4f_add_vec4f(vec4f_t first, vec4f_t second){
    vec4f_t res;
    for (int i = 0; i < 4; ++i) {
        ((float*)&res)[i] = ((float*)&first)[i] + ((float*)&second)[i];
    }
    return res;
}
vec4f_t vec4f_subtract_vec4f(vec4f_t first, vec4f_t second){
    vec4f_t res;
    for (int i = 0; i < 4; ++i) {
        ((float*)&res)[i] = ((float*)&first)[i] - ((float*)&second)[i];
    }
    return res;
}
vec4f_t vec4f_multiply_vec4f(vec4f_t first, vec4f_t second){
    vec4f_t res;
    for (int i = 0; i < 4; ++i) {
        ((float*)&res)[i] = ((float*)&first)[i] * ((float*)&second)[i];
    }
    return res;
}
vec4f_t vec4f_divide_vec4f(vec4f_t first, vec4f_t second){
    vec4f_t res;
    for (int i = 0; i < 4; ++i) {
        ((float*)&res)[i] = ((float*)&first)[i] / ((float*)&second)[i];
    }
    return res;
}
vec4f_t vec4f_add_float(vec4f_t vec, float value){
    vec4f_t res;
    for (int i = 0; i < 4; ++i) {
        ((float*)&res)[i] = ((float*)&vec)[i] + value;
    }
    return res;
}
vec4f_t vec4f_subtract_float(vec4f_t vec, float value){
    vec4f_t res;
    for (int i = 0; i < 4; ++i) {
        ((float*)&res)[i] = ((float*)&vec)[i] - value;
    }
    return res;
}
vec4f_t vec4f_multiply_float(vec4f_t vec, float value){
    vec4f_t res;
    for (int i = 0; i < 4; ++i) {
        ((float*)&res)[i] = ((float*)&vec)[i] * value;
    }
    return res;
}
vec4f_t vec4f_divide_float(vec4f_t vec, float value){
    vec4f_t res;
    for (int i = 0; i < 4; ++i) {
        ((float*)&res)[i] = ((float*)&vec)[i] / value;
    }
    return res;
}


/**** OPERATIONS  vec3f ****/
vec3f_t vec3f_multiply_mat3(vec3f_t vec, mat3_t mat){
    return VEC3F();
}
float   vec3f_dot_vec3f(vec3f_t first, vec3f_t second){
    return (first.x * second.x + first.y * second.y + first.z * second.z);
}
vec3f_t vec3f_cross_vec3f(vec3f_t first, vec3f_t second){
    vec3f_t res;

    res.x = first.y * second.z - first.z * second.y;
    res.y = first.z * second.x - first.x * second.z;
    res.z = first.x * second.y - first.y * second.x;

    return res;
}
float   vec3f_length(vec3f_t vec){
    return sqrtf(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}
vec3f_t vec3f_normalize_copy(vec3f_t vec){
    float l = vec3f_length(vec);
    if(l > 0)
        return vec3f(vec.x/l,vec.y/l,vec.z/l);
    else
        return vec3f(0,0,0);
}
void    vec3f_normalize(vec3f_t* vec){}
vec3f_t vec3f_add_vec3f(vec3f_t first, vec3f_t second){
    vec3f_t res;
    for (int i = 0; i < 3; ++i) {
        ((float*)&res)[i] = ((float*)&first)[i] + ((float*)&second)[i];
    }
    return res;
}
vec3f_t vec3f_subtract_vec3f(vec3f_t first, vec3f_t second){
    vec3f_t res;
    for (int i = 0; i < 3; ++i) {
        ((float*)&res)[i] = ((float*)&first)[i] - ((float*)&second)[i];
    }
    return res;
}
vec3f_t vec3f_multiply_vec3f(vec3f_t first, vec3f_t second){
    vec3f_t res;
    for (int i = 0; i < 3; ++i) {
        ((float*)&res)[i] = ((float*)&first)[i] * ((float*)&second)[i];
    }
    return res;
}
vec3f_t vec3f_divide_vec3f(vec3f_t first, vec3f_t second){
    vec3f_t res;
    for (int i = 0; i < 3; ++i) {
        ((float*)&res)[i] = ((float*)&first)[i] / ((float*)&second)[i];
    }
    return res;
}
vec3f_t vec3f_add_float(vec3f_t vec, float value){
    vec3f_t res;
    for (int i = 0; i < 3; ++i) {
        ((float*)&res)[i] = ((float*)&vec)[i] + value;
    }
    return res;
}
vec3f_t vec3f_subtract_float(vec3f_t vec, float value){
    vec3f_t res;
    for (int i = 0; i < 3; ++i) {
        ((float*)&res)[i] = ((float*)&vec)[i]  - value;
    }
    return res;
}
vec3f_t vec3f_multiply_float(vec3f_t vec, float value){
    vec3f_t res;
    for (int i = 0; i < 3; ++i) {
        ((float*)&res)[i] = ((float*)&vec)[i] * value;
    }
    return res;
}
vec3f_t vec3f_divide_float(vec3f_t vec, float value){
    vec3f_t res;
    for (int i = 0; i < 3; ++i) {
        ((float*)&res)[i] = ((float*)&vec)[i] / value;
    }
    return res;
}



/**** OPERATIONS  vec2f ****/
float   vec2f_dot_vec2f(vec2f_t first, vec2f_t second){
    return (first.x * second.x + first.y * second.y);
}
float   vec2f_length(vec2f_t vec){
    return sqrtf(vec.x*vec.x + vec.y*vec.y);
}
vec2f_t vec2f_normalize_copy(vec2f_t vec){
    float l = vec2f_length(vec);
    if(l > 0)
        return vec2f(vec.x/l,vec.y/l);
    else
        return vec2f(0,0);
}
void    vec2f_normalize(vec2f_t* vec){}
vec2f_t vec2f_add_vec2f(vec2f_t first, vec2f_t second){
    vec2f_t res;
    for (int i = 0; i < 2; ++i) {
        ((float*)&res)[i] = ((float*)&first)[i] + ((float*)&second)[i];
    }
    return res;
}
vec2f_t vec2f_subtract_vec2f(vec2f_t first, vec2f_t second){
    vec2f_t res;
    for (int i = 0; i < 2; ++i) {
        ((float*)&res)[i] = ((float*)&first)[i] - ((float*)&second)[i];
    }
    return res;
}
vec2f_t vec2f_multiply_vec2f(vec2f_t first, vec2f_t second){
    vec2f_t res;
    for (int i = 0; i < 2; ++i) {
        ((float*)&res)[i] = ((float*)&first)[i] * ((float*)&second)[i];
    }
    return res;
}
vec2f_t vec2f_divide_vec2f(vec2f_t first, vec2f_t second){
    vec2f_t res;
    for (int i = 0; i < 2; ++i) {
        ((float*)&res)[i] = ((float*)&first)[i] / ((float*)&second)[i];
    }
    return res;
}
vec2f_t vec2f_add_float(vec2f_t vec, float value){
    vec2f_t res;
    for (int i = 0; i < 2; ++i) {
        ((float*)&res)[i] = ((float*)&vec)[i] + value;
    }
    return res;
}
vec2f_t vec2f_subtract_float(vec2f_t vec, float value){
    vec2f_t res;
    for (int i = 0; i < 2; ++i) {
        ((float*)&res)[i] = ((float*)&vec)[i] - value;
    }
    return res;
}
vec2f_t vec2f_multiply_float(vec2f_t vec, float value){
    vec2f_t res;
    for (int i = 0; i < 2; ++i) {
        ((float*)&res)[i] = ((float*)&vec)[i] * value;
    }
    return res;
}
vec2f_t vec2f_divide_float(vec2f_t vec, float value){
    vec2f_t res;
    for (int i = 0; i < 2; ++i) {
        ((float*)&res)[i] = ((float*)&vec)[i] / value;
    }
    return res;
}



/**** OPERATIONS mat4 ****/
mat4_t mat4_multiply_mat4(mat4_t first, mat4_t second){
    mat4_t result;

    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            float sum = 0;
            for(int k = 0; k < 4; k++) {
                sum += first.elements[k][j] * second.elements[i][k];
            }
            result.elements[i][j] = sum;
        }
    }
    return result;
}
mat4_t mat4_invert_affine(mat4_t mat){
    mat4_t result;
     // Link : http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
    //TODO: finish
    float det = mat.m00 * mat.m11 * mat.m22 * mat.m33;
    det += mat.m00 * mat.m12 * mat.m23 * mat.m31;
    det += mat.m00 * mat.m13 * mat.m21 * mat.m32;

    det += mat.m01 * mat.m10 * mat.m23 * mat.m32;
    det += mat.m01 * mat.m12 * mat.m20 * mat.m33;
    det += mat.m01 * mat.m13 * mat.m22 * mat.m30;

    det += mat.m02 * mat.m10 * mat.m21 * mat.m33;
    det += mat.m02 * mat.m11 * mat.m23 * mat.m30;
    det += mat.m02 * mat.m13 * mat.m20 * mat.m31;

    return result;
}
mat4_t mat4_transpose(mat4_t matrix){
    return mat4(
            matrix.m00, matrix.m01, matrix.m02, matrix.m03,
            matrix.m10, matrix.m11, matrix.m12, matrix.m13,
            matrix.m20, matrix.m21, matrix.m22, matrix.m23,
            matrix.m30, matrix.m31, matrix.m32, matrix.m33
    );
}
mat4_t mat4_perspective_deg(float fovy, float near, float far, float aspect){
    return mat4_perspective_rad(DEG_TO_RAD(fovy),near,far,aspect);
}
mat4_t mat4_perspective_rad(float fovy, float near, float far, float aspect){
   // float fovy_in_rad = vertical_field_of_view_in_deg / 180 * M_PI;
    float f = 1.0f / tanf(fovy / 2.0f);
    float ar = aspect;
    float nd = near, fd = far;

    return mat4(
            f / ar,           0,                0,                0,
            0,                f,                0,                0,
            0,                0,               (fd+nd)/(nd-fd),  (2*fd*nd)/(nd-fd),
            0,                0,               -1,                0
    );
}
mat4_t mat4_look_at(vec3f_t from, vec3f_t to, vec3f_t up){
    vec3f_t z = vec3f_multiply_float(vec3f_normalize_copy(vec3f_subtract_vec3f(to, from)), -1);
    vec3f_t x = vec3f_normalize_copy(vec3f_cross_vec3f(up, z));
    vec3f_t y = vec3f_cross_vec3f(z, x);

    return mat4(
            x.x, x.y, x.z, -vec3f_dot_vec3f(from, x),
            y.x, y.y, y.z, -vec3f_dot_vec3f(from, y),
            z.x, z.y, z.z, -vec3f_dot_vec3f(from, z),
            0,   0,   0,    1
    );
}
mat4_t mat4_orthogonal(float left, float right, float top, float bottom, float back, float front){
    float l = left, r = right, b = bottom, t = top, n = front, f = back;
    float tx = -(r + l) / (r - l);
    float ty = -(t + b) / (t - b);
    float tz = -(f + n) / (f - n);
    return mat4(
            2 / (r - l),  0,            0,            tx,
            0,            2 / (t - b),  0,            ty,
            0,            0,            2 / (f - n),  tz,
            0,            0,            0,            1
    );
}
mat4_t mat4_rotation_deg(float angle,vec3f_t axis){
    return mat4_rotation_rad(DEG_TO_RAD(angle),axis);
}
mat4_t mat4_rotation_rad(float angle,vec3f_t axis){
    vec3f_t normalized_axis = vec3f_normalize_copy(axis);
    float x = normalized_axis.x, y = normalized_axis.y, z = normalized_axis.z;
    float c = cosf(angle), s = sinf(angle);

    return mat4(
            c + x*x*(1-c),            x*y*(1-c) - z*s,      x*z*(1-c) + y*s,  0,
            y*x*(1-c) + z*s,  c + y*y*(1-c),            y*z*(1-c) - x*s,  0,
            z*x*(1-c) - y*s,      z*y*(1-c) + x*s,  c + z*z*(1-c),        0,
            0,                        0,                    0,            1
    );;
}
mat4_t mat4_translation(vec3f_t translation) {
    return mat4(
            1,0,0,translation.x,
            0,1,0,translation.y,
            0,0,1,translation.z,
            0,0,0,1
    );
}
mat4_t mat4_scaling(vec3f_t scl){
    return mat4(
            scl.x,0,0,0,
            0,scl.y,0,0,
            0,0,scl.z,0,
            0,0,0,1
    );
}

/**** handy functions ****/
int    m3d_in_range(float value, float min, float max){
    if(value >= min && value <= max){
        return 1;
    }else{
        return 0;
    }
}
float  m3d_map(float value, float in_min, float in_max, float out_min, float out_max) {
    return (value - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}
float  m3d_constrain(float value, float min, float max){
    if(m3d_in_range(value,min,max)){
        return value;
    }
    if(value < min){
        return min;
    }
    if(value > max){
        return max;
    }
}
float  m3d_interpolate_linear(float start, float end, float factor){
    return (1 - factor) * start + factor * end;
}

/**** Debug ****/
void m3d_print_vec2f_raw(vec2f_t vec, const char* name){
    printf("[vec2f, %s]: (%f,%f)\n",name,vec.x,vec.y);
}
void m3d_print_vec3f_raw(vec3f_t vec, const char* name){
    printf("[vec3f, %s]: (%f,%f,%f)\n",name,vec.x,vec.y,vec.z);
}
void m3d_print_vec4f_raw(vec4f_t vec, const char* name){
    printf("[vec4f, %s]: (%f,%f,%f,%f)\n",name,vec.x,vec.y,vec.z,vec.w);
}
void m3d_print_mat4_raw(mat4_t mat, const char* name){
    printf("[mat4, %s \n]:",name);
    for (int i = 0; i < 4; ++i) {
        printf("\t(%f,%f,%f,%f)\n",mat.elements[i][0],mat.elements[i][1],mat.elements[i][2],mat.elements[i][3]);
    }
}
void m3d_print_mat3_raw(mat3_t mat, const char* name){
    printf("[mat3, %s]: \n",name);
    for (int i = 0; i < 3; ++i) {
        printf("\t(%f,%f,%f)\n",mat.elements[i][0],mat.elements[i][1],mat.elements[i][2]);
    }
}
