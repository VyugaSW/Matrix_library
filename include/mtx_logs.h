#pragma once

#include <stdio.h>

/**
 * @def MTX_DEBUG
 * @brief Debug mode switch (1 = enabled, 0 = disabled)
 * @details When enabled, all logging operations will be active. 
 * When disabled, all logging macros become no-ops for zero runtime overhead.
 */
#define MTX_DEBUG 1

/**
 * @def MTX_LOG_FILE
 * @brief Path to the log file
 * @note The directory must exist before logging starts
 */
#define MTX_LOG_FILE "data/mtx_log.txt"

#if MTX_DEBUG
    /**
     * @def MTX_LOG(msg)
     * @brief Logs an informational message
     * @param msg The message string to log
     * @details Writes to both log file and stdout. 
     * Creates/append to the file specified in MTX_LOG_FILE.
     */
    #define MTX_LOG(msg) do { \
        FILE *f = fopen(MTX_LOG_FILE, "a"); \
        if (f) { \
            fprintf(f, "[MTX] %s\n", msg); \
            fclose(f); \
        } \
    } while (0)

    /**
     * @def MTX_LOG_ERROR(msg)
     * @brief Logs an error message
     * @param msg The error message to log
     * @details Writes to both log file and stderr.
     * Creates/append to the file specified in MTX_LOG_FILE.
     */
    #define MTX_LOG_ERROR(msg) do { \
        FILE *f = fopen(MTX_LOG_FILE, "a"); \
        if (f) { \
            fprintf(f, "[MTX_ERROR] %s\n", msg); \
            fclose(f); \
        } \
        fprintf(stderr, "[MTX_ERROR] %s\n", msg); \
    } while (0)
#else
    /**
     * @def MTX_LOG(msg)
     * @brief Empty macro when debugging is disabled
     */
    #define MTX_LOG(msg) 
    
    /**
     * @def MTX_LOG_ERROR(msg)
     * @brief Empty macro when debugging is disabled
     */
    #define MTX_LOG_ERROR(msg)
#endif