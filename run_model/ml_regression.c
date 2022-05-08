/******************************************************************************
* File Name:   ml_local_regression.c
*
* Description: This file contains the implementation of local validation of the
*   machine learning model.
*
* Related Document: See README.md
*
*
*******************************************************************************
* Copyright 2021, Cypress Semiconductor Corporation (an Infineon company) or
* an affiliate of Cypress Semiconductor Corporation.  All rights reserved.
*
* This software, including source code, documentation and related
* materials ("Software") is owned by Cypress Semiconductor Corporation
* or one of its affiliates ("Cypress") and is protected by and subject to
* worldwide patent protection (United States and foreign),
* United States copyright laws and international treaty provisions.
* Therefore, you may use this Software only as provided in the license
* agreement accompanying the software package from which you
* obtained this Software ("EULA").
* If no EULA applies, Cypress hereby grants you a personal, non-exclusive,
* non-transferable license to copy, modify, and compile the Software
* source code solely for use in connection with Cypress's
* integrated circuit products.  Any reproduction, modification, translation,
* compilation, or representation of this Software except as specified
* above is prohibited without the express written permission of Cypress.
*
* Disclaimer: THIS SOFTWARE IS PROVIDED AS-IS, WITH NO WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, NONINFRINGEMENT, IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. Cypress
* reserves the right to make changes to the Software without notice. Cypress
* does not assume any liability arising out of the application or use of the
* Software or any product or circuit described in the Software. Cypress does
* not authorize its products for use in any products where a malfunction or
* failure of the Cypress product may reasonably be expected to result in
* significant property damage, injury or death ("High Risk Product"). By
* including Cypress's product in a High Risk Product, the manufacturer
* of such system or application assumes all risk of such use and in doing
* so agrees to indemnify Cypress against all liability.
*******************************************************************************/
#include "ml_local_regression.h"
#include "mtb_ml_utils.h"

#include "cyhal.h"

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

/* Include regression files */
#include MTB_ML_INCLUDE_MODEL_X_DATA_FILE(MODEL_NAME)
#include MTB_ML_INCLUDE_MODEL_Y_DATA_FILE(MODEL_NAME)

/*******************************************************************************
* Typedefs
*******************************************************************************/

/*******************************************************************************
* Constants
*******************************************************************************/
#define SUCCESS_RATE       98.0f

/*******************************************************************************
* Global Variables
*******************************************************************************/
/* NN Model Object */
static mtb_ml_model_t *model_obj;

/* Output/result buffers for the inference engine */
static MTB_ML_DATA_T *result_buffer;

/* Model Output Size */
static uint32_t model_output_size;

/*******************************************************************************
* Local Functions
*******************************************************************************/

/*******************************************************************************
* Function Name: nn_profiler_init
********************************************************************************
* Summary:
*   Initialize the Neural Network based on the given model and setup to start
*   regression of the model and profiling configuration.
*
* Parameters:
*   profile_cfg: profiling configuration
*   model_bin: pointer to the model data
*
* Return:
*   The status of the initialization.
*******************************************************************************/
cy_rslt_t ml_local_regression_init(mtb_ml_profile_config_t profile_cfg,
                                   mtb_ml_model_bin_t *model_bin)
{
    cy_rslt_t result;

    /* Initialize the neural network */
    result = mtb_ml_model_init(model_bin,
                               NULL,
                               &model_obj);

    if (result != CY_RSLT_SUCCESS)
    {
        printf("MTB ML initialization failure: %lu\r\n", (unsigned long) result);
        return result;
    }

    mtb_ml_model_profile_config(model_obj, profile_cfg);

    model_output_size = mtb_ml_model_get_output_size(model_obj);
    printf("%d \r\n", sizeof(MTB_ML_DATA_T));

    /* Allocate memory for the output buffers */
    result_buffer = (MTB_ML_DATA_T *) malloc(16);
    // result_buffer = (MTB_ML_DATA_T *) malloc(model_output_size * sizeof(MTB_ML_DATA_T));
    if (result_buffer == NULL)
    {
        return MTB_ML_RESULT_ALLOC_ERR;
    };

    /* Print information about the model */
    mtb_ml_model_info(model_obj, model_bin);
    
    return CY_RSLT_SUCCESS;
}


void ml_task(float spectogram[2816]) {
	int16_t spectogram_bits[2816];
	mtb_ml_utils_convert_flt_to_int16(spectogram, spectogram_bits, 2816, 3.13);

	MTB_ML_DATA_T *in_buff = spectogram_bits;
	cy_rslt_t msg = mtb_ml_model_set_input_q_fraction_bits(model_obj, 3.13);

	mtb_ml_model_run(model_obj, in_buff, result_buffer);
	int i;
	for (i = 0; i < 8; i++) {
		printf("%d: %d\r\n", i, result_buffer[i]);
	}

	i = mtb_ml_utils_find_max(result_buffer, model_output_size);
	printf("max: %d\r\n", i);

}