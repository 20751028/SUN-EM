// ****************************************
// Optix Raytracer Interface for Matlab
//
// mexRaytracer.cpp
//
// CD Stewart-Burger 16/01/2020
// 
//
// ****************************************

// Global defines
#define MAX_ARG_LEN 256
#define MAX_MSG_LEN 256

// Stop windows.h from defining min and max as required by optix_math.h
#define NOMINMAX

// Default includes
#include "mex.h"
#include "optix_world.h"
#include "optix_prime.h"
#include <optixu/optixu_math_namespace.h>
#include <string.h>

typedef optix::float3 float3;
typedef optix::int3 int3;

struct Ray
{
    float3 origin;
    float3 direction;
};

struct Hit
{
    float t;
    int triId;
};

// Platform specific includes
#ifdef WIN32
#include <process.h>
#else
#include <pthread.h>
#endif

// Global variables
RTPcontext g_context = NULL;
RTPresult g_code = RTP_SUCCESS;
RTPmodel g_model;
RTPbuffertype g_bufferType = RTP_BUFFER_TYPE_HOST;
RTPbufferdesc raysBD;
RTPbufferdesc hitsBD;
Hit *g_hits = NULL;
Ray *g_rays = NULL;
int g_numthreads = 0;
int3 *g_indexdata = NULL;
float3 *g_vertexdata = NULL;
unsigned int g_numrays = 0;
char g_time_str[MAX_ARG_LEN + 1];

bool bInitRequired = true;
bool bRunNonBlocking = false;
bool bResultsAvail = false;
bool bResultsRead = false;
bool bVerbose = false;

#ifdef WIN32
HANDLE g_thread;
#else
pthread_t g_thread;
#endif

#ifdef WIN32
void RunTraversal(void *)
#else
void *RunTraversal(void *)
#endif
{

    RTPquery query;
    g_code = RTP_ERROR_UNKNOWN;
    if (rtpQueryCreate(g_model, RTP_QUERY_TYPE_CLOSEST, &query) != RTP_SUCCESS)
    {
        mexWarnMsgTxt("rtpQueryCreate failed\n");
    }
    else if (rtpQuerySetRays(query, raysBD) != RTP_SUCCESS)
    {
        mexWarnMsgTxt("rtpQuerySetRays failed\n");
    }
    else if (rtpQuerySetHits(query, hitsBD) != RTP_SUCCESS)
    {
        mexWarnMsgTxt("rtpQuerySetHits failed\n");
    }
    else
    {
        rtpModelFinish(g_model);
        if (rtpQueryExecute(query, 0 /* hints */) != RTP_SUCCESS)
        {
            mexWarnMsgTxt("rtpQueryExecute failed\n");
        }
        else
            g_code = RTP_SUCCESS;
        rtpQueryDestroy(query);
    }
}

void WaitForThreadToFinish(void)
{
    if (bVerbose)
        mexPrintf("Waiting for traversal thread to finish\n");

        // Wait for traversal thread to finish
#ifdef WIN32
    WaitForSingleObject(g_thread, INFINITE);
#else
    pthread_join(g_thread, NULL);
#endif
    bRunNonBlocking = false;
    bResultsAvail = true;
    bResultsRead = false;
    if (g_code != RTP_SUCCESS)
    {
        bResultsAvail = false;
        mexWarnMsgTxt("Non-blocking run of rtuTraversalTraverse failed\n");
    }
}

void Initialise(void)
{
    if (bVerbose)
        mexPrintf("Initializing\n");
    // Create traversal API & context
    if (rtpContextCreate(RTP_CONTEXT_TYPE_CUDA, &g_context) != RTP_SUCCESS)
    {
        if (rtpContextCreate(RTP_CONTEXT_TYPE_CPU, &g_context) != RTP_SUCCESS)
        {
            mexErrMsgTxt("rtpContextCreate failed\n");
        }
        else
        {
            if (g_numthreads == 0)
            {
                g_numthreads = 1;
            }
            if (rtpContextSetCpuThreads(g_context, g_numthreads) != RTP_SUCCESS)
            {
                mexErrMsgTxt("rtpContextSetCpuThreads failed\n");
            }
            else
                mexPrintf("Using CPU (GPU context setup failed)\n");
        }
        g_bufferType = RTP_BUFFER_TYPE_HOST;
        //*
    }
    else
    {
        if (bVerbose)
            mexPrintf("Using GPU\n");
        //g_bufferType = RTP_BUFFER_TYPE_CUDA_LINEAR;//change to CUDA type
    } //*/
}

static void CleanUp(void)
{
    if (bVerbose)
        mexPrintf("Cleaning up\n");

    // Wait for previous run to finish
    if (bRunNonBlocking)
    {
        WaitForThreadToFinish();
    }

    // Destroy traversal API & context
    g_code = rtpContextDestroy(g_context);

    // Make sure memory gets freed
    if (g_vertexdata)
    {
       // mexPrintf("Freeing vertices\n");
        mxFree(g_vertexdata);
    }
    if (g_indexdata)
    {
       // mexPrintf("Freeing indices\n");
        mxFree(g_indexdata);
    }
    if (g_rays)
    {
       // mexPrintf("Freeing rays\n");
        mxFree(g_rays);
    }
    if (g_hits)
    {
        //mexPrintf("Freeing hits\n");
        mxFree(g_hits);
    }

    // Now its ok to throw error, if any
    if (g_code != RTP_SUCCESS)
    {
        mexErrMsgTxt("rtpContextDestroy failed\n");
    }
}

template <class T_src>
void transposeFloat3(float3 *p_dst, T_src *p_src,
                     unsigned int R_src)
{
    unsigned int r;
    for (r = 0; r < R_src; r++)
    {
        p_dst[r].x = (float)p_src[r];
        p_dst[r].y = (float)p_src[r + R_src];
        p_dst[r].z = (float)p_src[r + 2 * R_src];
    }
}

template <class T_src>
void transposeInt3(int3 *p_dst, T_src *p_src,
                   unsigned int R_src)
{
    unsigned int r;
    for (r = 0; r < R_src; r++)
    {
        p_dst[r].x = (int)p_src[r];
        p_dst[r].y = (int)p_src[r + R_src];
        p_dst[r].z = (int)p_src[r + 2 * R_src];
    }
}

template <class T_src>
void transposeRay(Ray *p_dst, T_src *p_src,
                  unsigned int R_src)
{
    unsigned int r;
    for (r = 0; r < R_src; r++)
    {
        p_dst[r].origin.x = (float)p_src[r];
        p_dst[r].origin.y = (float)p_src[r + R_src];
        p_dst[r].origin.z = (float)p_src[r + 2 * R_src];
        p_dst[r].direction.x = (float)p_src[r + 3 * R_src];
        p_dst[r].direction.y = (float)p_src[r + 4 * R_src];
        p_dst[r].direction.z = (float)p_src[r + 5 * R_src];
    }
}

// Repack data from (x1 x2 x3 ... y1 y2 y3 ... z1 z2 z3 ...)
// to (x1 y1 z1 x2 y2 z2 ...) and convert types as required
template <class T_dst, class T_src>
void Transpose(T_dst *p_dst, T_src *p_src,
               unsigned int R_src, unsigned int C_src)
{
    unsigned int r, c;
    for (c = 0; c < C_src; c++)
    {
        for (r = 0; r < R_src; r++)
        {
            *p_dst++ = (T_dst)p_src[(r * C_src) + c];
        }
    }
}

// Convert data type from a source to destination block
template <class T_dst, class T_src>
void Convert(T_dst *p_dst, T_src *p_src, unsigned int N_src)
{
    unsigned int n;
    for (n = 0; n < N_src; n++)
    {
        *p_dst++ = (T_dst)p_src[n];
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char arg_str[MAX_ARG_LEN + 1];

    // Check inputs
    if (nrhs < 1)
        mexErrMsgTxt("No arguments\n");
    if (!mxIsChar(prhs[0]))
        mexErrMsgTxt("First argument must be a command string\n");
    mxGetString(prhs[0], arg_str, MAX_ARG_LEN);

    // Wait for previous run to finish
    if (bRunNonBlocking)
    {
        WaitForThreadToFinish();
    }

    // Handle CPU thread argument first (in case of CPU only init)
    if (strcmp(arg_str, "cputhreads") == 0)
    {
        // Check number of CPU threads argument
        if (nrhs < 2)
            mexErrMsgTxt("Number of CPU threads not specified\n");
        if (!mxIsNumeric(prhs[1]))
        {
            mexErrMsgTxt("Number of CPU threads must be an integer\n");
        }
        int numthreads = (int)(*mxGetPr(prhs[1]));
        if (numthreads > 0)
        {
            g_numthreads = numthreads;
        }
        else
        {
            mexErrMsgTxt("Invalid number of CPU threads\n");
        }

        if (!bInitRequired)
        {
            // Set number of CPU threads to use when GPU is not available
            // and ray traversal reverts to CPU implementation
            if (rtpContextSetCpuThreads(g_context, g_numthreads) != RTP_SUCCESS)
            {
                mexErrMsgTxt("rtpContextSetCpuThreads failed\n");
            }
        }
    }

    // Handle initialisation
    if (bInitRequired)
    {
        mexAtExit(CleanUp);
        Initialise();
        bInitRequired = 0;
    }

    // Process arguments
    if (strcmp(arg_str, "cputhreads") == 0)
    {
        // Detect argument as valid
    }
    else if (strcmp(arg_str, "setmesh") == 0)
    {
        // Check triangle mesh inputs
        if (nrhs < 2)
            mexErrMsgTxt("No vertex data specified\n");
        if (!(mxIsSingle(prhs[1]) || mxIsDouble(prhs[1])))
        {
            mexErrMsgTxt("Vertex data must be a single or double matrix\n");
        }
        if (mxGetNumberOfDimensions(prhs[1]) != 2)
        {
            mexErrMsgTxt("Vertex data must be a 2D matrix\n");
        }
        if ((mxGetM(prhs[1]) != 3) && (mxGetN(prhs[1]) != 3))
        {
            mexErrMsgTxt("Vertex data must be (x y z)*n_vertices\n");
        }
        unsigned int n_vertices = (unsigned int)mxGetNumberOfElements(prhs[1]) / 3;
        if (n_vertices == 0)
        {
            mexErrMsgTxt("Vertex data empty\n");
        }

        if (nrhs < 3)
            mexErrMsgTxt("No index data specified\n");
        if (!(mxIsUint32(prhs[2]) || mxIsDouble(prhs[2])))
        {
            mexErrMsgTxt("Index data must be a uint32 or double index matrix\n");
        }
        if (mxGetNumberOfDimensions(prhs[2]) != 2)
        {
            mexErrMsgTxt("Index data must be a 2D matrix\n");
        }
        if ((mxGetM(prhs[2]) != 3) && (mxGetN(prhs[2]) != 3))
        {
            mexErrMsgTxt("Index data must be (i1 i2 i3)*n_indices\n");
        }
        unsigned int n_indices = (unsigned int)mxGetNumberOfElements(prhs[2]) / 3;
        if (n_indices == 0)
        {
            mexErrMsgTxt("Index data empty\n");
        }

        // Compare the incoming geometry 'last modified' to the existing geometry
        if (nrhs < 4)
            mexErrMsgTxt("No Mesh time\n");
        if (!mxIsChar(prhs[3]))
            mexErrMsgTxt("Mesh time must be a string\n");
        char time_str[MAX_ARG_LEN + 1];
        mxGetString(prhs[3], time_str, MAX_ARG_LEN);
        if (strcmp(g_time_str, time_str) == 0)
        {
            return;
        }
        strcpy(g_time_str, time_str);

        // Copy vertex data
        if (g_vertexdata)
        {
            mxFree(g_vertexdata);
        }
        g_vertexdata = (float3 *)mxMalloc(n_vertices * sizeof(float3));
        if (g_vertexdata)
        {
            mexMakeMemoryPersistent(g_vertexdata);
        }
        else
        {
            mexErrMsgTxt("Unable to allocate memory for vertex data\n");
        }
        {
            // Tranpose and/or convert vertex data as required
            float3 *p_dst = g_vertexdata;
            if (mxIsSingle(prhs[1]))
            {
                float *p_src = (float *)mxGetData(prhs[1]);
                transposeFloat3<float>(p_dst, p_src, n_vertices);
            }
            else
            {
                double *p_src = mxGetPr(prhs[1]);
                transposeFloat3<double>(p_dst, p_src, n_vertices);
            }
        }

        // Copy index data
        if (g_indexdata)
        {
            mxFree(g_indexdata);
        }
        g_indexdata = (int3 *)mxMalloc(n_indices * sizeof(int3));
        if (g_indexdata)
        {
            mexMakeMemoryPersistent(g_indexdata);
        }
        else
        {
            mexErrMsgTxt("Unable to allocate memory for index data\n");
        }
        {
            // Tranpose and/or convert index data as required
            int3 *p_dst = g_indexdata;
            if (mxIsUint32(prhs[2]))
            {
                unsigned int *p_src = (unsigned int *)mxGetData(prhs[2]);
                transposeInt3<unsigned int>(p_dst, p_src, n_indices);
            }
            else
            {
                // Convert and interleave data
                double *p_src = mxGetPr(prhs[2]);
                transposeInt3<double>(p_dst, p_src, n_indices);
            }
        }

        RTPbufferdesc verticesBD;
        if (rtpBufferDescCreate(g_context, RTP_BUFFER_FORMAT_VERTEX_FLOAT3, RTP_BUFFER_TYPE_HOST, g_vertexdata, &verticesBD) != RTP_SUCCESS)
        {
            mexErrMsgTxt("rtpBufferDescCreate failed for vertexdata\n");
        }
        else if (rtpBufferDescSetRange(verticesBD, 0, n_vertices) != RTP_SUCCESS)
        {
            mexErrMsgTxt("rtpBufferDescSetRange failed for vertexdata\n");
        }

        RTPbufferdesc indicesBD;
        if (rtpBufferDescCreate(g_context, RTP_BUFFER_FORMAT_INDICES_INT3, RTP_BUFFER_TYPE_HOST, g_indexdata, &indicesBD) != RTP_SUCCESS)
        {
            mexErrMsgTxt("rtpBufferDescCreate failed\n");
        }
        else if (rtpBufferDescSetRange(indicesBD, 0, n_indices) != RTP_SUCCESS)
        {
            mexErrMsgTxt("rtpBufferDescSetRange failed for indexdata\n");
        }
        // Set traingle mesh data
        if (g_model)
            rtpModelFinish(g_model);
        if (rtpModelCreate(g_context, &g_model) != RTP_SUCCESS)
        {
            mexErrMsgTxt("rtpModelCreate failed\n");
        }
        if (rtpModelSetTriangles(g_model, indicesBD, verticesBD) != RTP_SUCCESS)
        {
            mexErrMsgTxt("rtpModelSetTriangles failed\n");
        }

        //clock_t start = clock();
        //Model update is curently done in blocking mode. Asynchronous mode is possible.
        if (rtpModelUpdate(g_model, RTP_MODEL_HINT_ASYNC) != RTP_SUCCESS)
        {
            mexErrMsgTxt("rtpModelUpdate failed\n");
        }

        //clock_t stop = clock();
        //mexPrintf("%lf", (double)(stop-start));
    }
    else if (strcmp(arg_str, "castrays") == 0)
    {
        // Check ray data inputs
        if (nrhs < 2)
            mexErrMsgTxt("Ray data missing\n");
        if (!(mxIsSingle(prhs[1]) || mxIsDouble(prhs[1])))
        {
            mexErrMsgTxt("Ray data must be a single or double matrix\n");
        }
        if (mxGetNumberOfDimensions(prhs[1]) != 2)
        {
            mexErrMsgTxt("Ray data must be a 2D matrix\n");
        }
        if ((mxGetM(prhs[1]) != 6) && (mxGetN(prhs[1]) != 6))
        {
            mexErrMsgTxt("Ray data must be (px py pz kx ky kz)*n_rays\n");
        }
        if (bResultsAvail && !bResultsRead)
        {
            mexWarnMsgTxt("Discarding previous results\n");
        }
        g_numrays = (unsigned int)mxGetNumberOfElements(prhs[1]) / 6;
        bResultsAvail = false;
        bResultsRead = false;
        if(g_numrays == 0) return;

        // Check blocking mode input
        if (nrhs == 3)
        {
            char mode_str[MAX_ARG_LEN + 1];

            if (!mxIsChar(prhs[2]))
                mexErrMsgTxt("Blocking mode must be specified with a string\n");
            mxGetString(prhs[2], mode_str, MAX_ARG_LEN);
            if (strcmp(mode_str, "noblock") == 0)
            {
                bRunNonBlocking = true;
            }
        }

        // Copy ray data
        {
            // Map ray buffer
            if (g_rays)
            {
                rtpBufferDescDestroy(raysBD);
                mxFree(g_rays);
            }
            g_rays = (Ray *)mxMalloc(g_numrays * sizeof(Ray));
            if (g_rays)
            {
                mexMakeMemoryPersistent(g_rays);
            }
            else
            {
                mexErrMsgTxt("Unable to allocate memory for rays\n");
            }

            // Tranpose and/or convert ray data as required
            if (mxIsSingle(prhs[1]))
            {
                float *p_src = (float *)mxGetData(prhs[1]);
                transposeRay<float>(g_rays, p_src, g_numrays);
            }
            else
            {
                double *p_src = mxGetPr(prhs[1]);

                transposeRay<double>(g_rays, p_src, g_numrays);
            }

            if (rtpBufferDescCreate(g_context, RTP_BUFFER_FORMAT_RAY_ORIGIN_DIRECTION, g_bufferType, g_rays, &raysBD) != RTP_SUCCESS)
            {
                mexErrMsgTxt("rtpBufferDescCreate failed for rays\n");
            }
            else if (rtpBufferDescSetRange(raysBD, 0, g_numrays) != RTP_SUCCESS)
            {
                mexErrMsgTxt("rtpBufferDescSetRange failed for rays\n");
            }
        }

        // Create buffer for returned hit descriptions
        if (g_hits)
        {
            rtpBufferDescDestroy(hitsBD);
            mxFree(g_hits);
        }
        g_hits = (Hit *)mxMalloc(g_numrays * sizeof(Hit));
        if (g_hits)
        {
            mexMakeMemoryPersistent(g_hits);
        }
        else
        {
            mexErrMsgTxt("Unable to allocate memory for hits\n");
        }
        if (rtpBufferDescCreate(g_context, RTP_BUFFER_FORMAT_HIT_T_TRIID, g_bufferType, g_hits, &hitsBD) != RTP_SUCCESS)
        {
            mexErrMsgTxt("rtpBufferDescCreate failed for hits\n");
        }
        else if (rtpBufferDescSetRange(hitsBD, 0, g_numrays) != RTP_SUCCESS)
        {
            mexErrMsgTxt("rtpBufferDescSetRange failed for vertexdata\n");
        }

        // Run traversal to cast rays
        if (bRunNonBlocking)
        {
            // Run traversal in separate thread (non-block)
#ifdef WIN32
            g_thread = (HANDLE)_beginthread(RunTraversal, 0, NULL);
#else
            pthread_create(&g_thread, NULL, &RunTraversal, NULL);
#endif
        }
        else
        {
            // Run traversal now and wait (block)
            RunTraversal(NULL);
            bResultsAvail = true;
            bResultsRead = false;
            if (g_code != RTP_SUCCESS)
            {
                bResultsAvail = false;
                mexErrMsgTxt("Traversal failed\n");
            }
        }
    }
    else if (strcmp(arg_str, "getresults") == 0)
    {
        // Check outputs
        if (nlhs < 1)
            mexErrMsgTxt("Primitive ID output missing\n");
        if (nlhs < 2)
            mexErrMsgTxt("Intersection distance output missing\n");

        // Check results
        if (!bResultsAvail){
            mexWarnMsgTxt("No results available\n");
        }
        plhs[0] = mxCreateNumericMatrix(g_numrays, 1, mxINT32_CLASS, mxREAL);
        int *p_prim_id = (int *)mxGetData(plhs[0]);

        plhs[1] = mxCreateNumericMatrix(g_numrays, 1, mxSINGLE_CLASS, mxREAL);
        float *p_dist = (float *)mxGetData(plhs[1]);
        /*
        The ray reported ray distance needs to be no greater than the actual ray distance
        to ensure that the reflected rays are not generated behind the triangle
        from which they are reflected. This requirement can be met by multiplying
        the ray distance by 1-(2^-23).
        */
        float m = 1.0f-pow(2, -23);
        for (unsigned int n = 0; n < g_numrays; n++)
        {
            *p_prim_id++ = g_hits[n].triId;
            *p_dist++ = g_hits[n].t*m;
        }
        bResultsRead = true;
    }
    else if (strcmp(arg_str, "getsurfnorm") == 0)
    {
        // Check output
        if (nlhs < 1)
            mexErrMsgTxt("Surface normal output missing\n");

        // Check results
        if (!bResultsAvail)
            mexErrMsgTxt("No results available\n");

        // Create output
        plhs[0] = mxCreateNumericMatrix(g_numrays, 3, mxSINGLE_CLASS, mxREAL);
        float *p_surfnorm = (float *)mxGetData(plhs[0]);

        float3 *norms = (float3 *)mxMalloc(g_numrays * sizeof(float3));
        if (norms)
        {
            //mexMakeMemoryPersistent(norms);
        }
        else
        {
            mexErrMsgTxt("Unable to allocate memory for norms\n");
        }

        for (size_t i = 0; i < g_numrays; i++)
        {
            if (g_hits[i].t < 0.0f)
            {
                //no hit
                norms[i].x = 0;
                norms[i].y = 0;
                norms[i].z = 0;
            }
            else
            {
                int3 tri = g_indexdata[g_hits[i].triId];
                float3 v0 = g_vertexdata[tri.x];
                float3 v1 = g_vertexdata[tri.y];
                float3 v2 = g_vertexdata[tri.z];
                float3 e0 = v1 - v0;
                float3 e1 = v2 - v0;
                norms[i] = normalize(cross(e0, e1));
            }
        }

        // Repack data from (nx1 ny1 nz1 nx2 ny2 nz2 ...)
        // to (nx1 nx2 nx3 ... ny1 ny2 ny3 ... nz1 nz2 nz3 ...)
        for (unsigned int n = 0; n < (g_numrays); n++)
        {
            *p_surfnorm++ = norms[n].x;
        }
        for (unsigned int n = 0; n < (g_numrays); n++)
        {
            *p_surfnorm++ = norms[n].y;
        }
        for (unsigned int n = 0; n < (g_numrays); n++)
        {
            *p_surfnorm++ = norms[n].z;
        }

        mxFree(norms);
    } //*/

    else if (strcmp(arg_str, "verbose") == 0)
    {
        // Set verbose option
        bVerbose = true;
    }
    else
    {
        char err_msg[MAX_MSG_LEN];
        sprintf(err_msg, "Unknown argument '%s'\n", arg_str);
        mexErrMsgTxt(err_msg);
    }
}

// EOF
