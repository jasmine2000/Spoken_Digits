/******************************************************************************
* File Name:   main.c
*
* Description: This is the source code for the Neural Network Profiler Example
*              for ModusToolbox.
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

#include "cy_pdl.h"
#include "cyhal.h"
#include "cybsp.h"
#include "cy_retarget_io.h"

#include "stdlib.h"

#include "elapsed_timer.h"
#include "mtb_ml_stream.h"
#include "ml_local_regression.h"

#include "spectogram.h"

#include MTB_ML_INCLUDE_MODEL_FILE(MODEL_NAME)

/*******************************************************************************
* Macros
********************************************************************************/
#define PROFILE_CONFIGURATION CY_ML_PROFILE_ENABLE_MODEL

/* Define how many samples in a frame */
#define FRAME_SIZE                  (512)
/* Noise threshold hysteresis */
#define THRESHOLD_HYSTERESIS        3u
/* Volume ratio for noise and print purposes */
#define VOLUME_RATIO                (4*FRAME_SIZE)
/* Desired sample rate. Typical values: 8/16/22.05/32/44.1/48kHz */
#define SAMPLE_RATE_HZ              8000u
/* Decimation Rate of the PDM/PCM block. Typical value is 64 */
#define DECIMATION_RATE             64u
/* Audio Subsystem Clock. Typical values depends on the desire sample rate:
- 8/16/48kHz    : 24.576 MHz
- 22.05/44.1kHz : 22.579 MHz */
#define AUDIO_SYS_CLOCK_HZ          24576000u
/* PDM/PCM Pins */
#define PDM_DATA                    P10_5
#define PDM_CLK                     P10_4

#define standby 0
#define recording 1
#define processing 2

/*******************************************************************************
* Function Prototypes
********************************************************************************/
void pdm_pcm_isr_handler(void *arg, cyhal_pdm_pcm_event_t event);
void clock_init(void);

/*******************************************************************************
* Global Variables
********************************************************************************/
/* Interrupt flags */
volatile bool pdm_pcm_flag = true;

/* Volume variables */
uint32_t volume = 0;
uint32_t noise_threshold = THRESHOLD_HYSTERESIS;

/* HAL Object */
cyhal_pdm_pcm_t pdm_pcm;
cyhal_clock_t   audio_clock;
cyhal_clock_t   pll_clock;

/* HAL Config */
const cyhal_pdm_pcm_cfg_t pdm_pcm_cfg =
{
    .sample_rate     = SAMPLE_RATE_HZ,
    .decimation_rate = DECIMATION_RATE,
    .mode            = CYHAL_PDM_PCM_MODE_STEREO,
    .word_length     = 16,  /* bits */
    .left_gain       = 0,   /* dB */
    .right_gain      = 0,   /* dB */
};

/* variables */
#define arr_length 5632
#define spect_length 2816


int recording_index;
int state;
int pressed;	// button

//MTB_ML_DATA_T *result_buffer;

float *spectogram;
double *voice_record;

/*******************************************************************************
* Function Name: main
********************************************************************************
* Summary:
* This is the main function for CM4 CPU. It does...
*    1. Initializes the BSP.
*    2. Prints welcome message
*    3. Initialize the regression unit - stream or local
*    4. Run the regression
*
* Parameters:
*  void
*
* Return:
*  int
*
*******************************************************************************/
int main(void)
{

    cy_rslt_t result;
    int16_t audio_frame[FRAME_SIZE] = {0};

    /* Initialize the device and board peripherals */
    result = cybsp_init() ;
    if (result != CY_RSLT_SUCCESS)
    {
        CY_ASSERT(0);
    }

    /* Enable global interrupts */
    __enable_irq();

    /* Init the clocks */
    clock_init();

    /* Initialize retarget-io to use the debug UART port */
	cy_retarget_io_init(CYBSP_DEBUG_UART_TX, CYBSP_DEBUG_UART_RX, CY_RETARGET_IO_BAUDRATE);
	printf("Initializing system...\r\n");

	/* Initialize ML model */
	mtb_ml_model_bin_t model_bin = {MTB_ML_MODEL_BIN_DATA(MODEL_NAME)};
	result = ml_local_regression_init(PROFILE_CONFIGURATION, &model_bin);
//	result_buffer = (MTB_ML_DATA_T *) malloc(16);
	printf("result: %lu", result);
	if (result == CY_RSLT_SUCCESS) printf("Successfully initialized model.");

	/* Initialize the User LED */
	cyhal_gpio_init(CYBSP_USER_LED, CYHAL_GPIO_DIR_OUTPUT, CYHAL_GPIO_DRIVE_STRONG, CYBSP_LED_STATE_OFF);

	/* Initialize the User Button */
	cyhal_gpio_init(CYBSP_USER_BTN, CYHAL_GPIO_DIR_INPUT, CYHAL_GPIO_DRIVE_PULLUP, CYBSP_BTN_OFF);

	/* Initialize the PDM/PCM block */
	cyhal_pdm_pcm_init(&pdm_pcm, PDM_DATA, PDM_CLK, &audio_clock, &pdm_pcm_cfg);
	cyhal_pdm_pcm_register_callback(&pdm_pcm, pdm_pcm_isr_handler, NULL);
	cyhal_pdm_pcm_enable_event(&pdm_pcm, CYHAL_PDM_PCM_ASYNC_COMPLETE, CYHAL_ISR_PRIORITY_DEFAULT, true);
	cyhal_pdm_pcm_start(&pdm_pcm);

	/* \x1b[2J\x1b[;H - ANSI ESC sequence for clear screen */
	printf("\x1b[2J\x1b[;H");

	pressed = 0;
	recording_index = 0;
	noise_threshold = 8;

	voice_record = (double *) malloc(arr_length);

	printf("\n\r");
    for (;;)
    {
    	/* Check if any microphone has data to process */
		if (pdm_pcm_flag)
		{
			/* Clear the PDM/PCM flag */
			pdm_pcm_flag = 0;

			pressed = !(cyhal_gpio_read(CYBSP_USER_BTN));

			/* Reset the volume */
			volume = 0;

			if (pressed || state == recording) {
				/* Calculate the volume by summing the absolute value of all the
				 * audio data from a frame */
				for (uint32_t index = 0; index < FRAME_SIZE; index++)
				{
					volume += abs(audio_frame[index]);
					voice_record[recording_index * FRAME_SIZE + index] = abs(audio_frame[index]);
				}

			}


			/* Turn ON the LED when the volume is higher than the threshold */
			if ((volume/VOLUME_RATIO) > noise_threshold || state == recording)
			{
				recording_index++;
				if (recording_index > 3) state = recording;

				cyhal_gpio_write((cyhal_gpio_t) CYBSP_USER_LED, CYBSP_LED_STATE_ON);
			}
			else
			{
				if (state == recording) {
					state = processing;
				} else {
					recording_index = 0;
				}

				cyhal_gpio_write((cyhal_gpio_t) CYBSP_USER_LED, CYBSP_LED_STATE_OFF);
			}

			/* Setup to read the next frame */
			cyhal_pdm_pcm_read_async(&pdm_pcm, audio_frame, FRAME_SIZE);
		}

		if (recording_index > 10) {
			state = processing;
		}

		if (state == processing) {
			cyhal_system_sleep();
			cyhal_pdm_pcm_stop(&pdm_pcm);

			int i;
			for (i = recording_index * FRAME_SIZE; i < arr_length; i++) {
				voice_record[i] = 0;
			}

			spectogram = (float *) malloc(spect_length);
			stft(voice_record, spectogram);

//			for (i = 0; i < 10; i++) {
//				printf("%f\r\n", spectogram[i]);
//			}

			free(voice_record);

			ml_task(spectogram);

			free(spectogram);

			recording_index = 0;
			state = standby;
			voice_record = (double *) malloc(arr_length);

			cyhal_pdm_pcm_start(&pdm_pcm);
			cyhal_pdm_pcm_read_async(&pdm_pcm, audio_frame, FRAME_SIZE);
		}

		cyhal_system_sleep();
    }
}

/*******************************************************************************
* Function Name: pdm_pcm_isr_handler
********************************************************************************
* Summary:
*  PDM/PCM ISR handler. Set a flag to be processed in the main loop.
*
* Parameters:
*  arg: not used
*  event: event that occurred
*
*******************************************************************************/
void pdm_pcm_isr_handler(void *arg, cyhal_pdm_pcm_event_t event)
{
    (void) arg;
    (void) event;

    pdm_pcm_flag = true;
}

/*******************************************************************************
* Function Name: clock_init
********************************************************************************
* Summary:
*  Initialize the clocks in the system.
*
*******************************************************************************/
void clock_init(void)
{
    /* Initialize the PLL */
    cyhal_clock_get(&pll_clock, &CYHAL_CLOCK_PLL[0]);
    cyhal_clock_init(&pll_clock);
    cyhal_clock_set_frequency(&pll_clock, AUDIO_SYS_CLOCK_HZ, NULL);
    cyhal_clock_set_enabled(&pll_clock, true, true);

    /* Initialize the audio subsystem clock (CLK_HF[1])
     * The CLK_HF[1] is the root clock for the I2S and PDM/PCM blocks */
    cyhal_clock_get(&audio_clock, &CYHAL_CLOCK_HF[1]);
    cyhal_clock_init(&audio_clock);

    /* Source the audio subsystem clock from PLL */
    cyhal_clock_set_source(&audio_clock, &pll_clock);
    cyhal_clock_set_enabled(&audio_clock, true, true);
}

/* [] END OF FILE */
