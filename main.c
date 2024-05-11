#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <wfdb/wfdb.h>
#include <wfdb/ecgmap.h>

typedef enum {false, true} bool;

#define WINDOWSIZE 20 
#define FS 360
#define BUFFSIZE 600 
#define MAX_LEADS 12

#define true 1
#define false 0

static WFDB_Siginfo s[MAX_LEADS];

unsigned stklen = 40000;
long nrSamps;
double *sig0;
double *sig4;

int m;
double dsum, dmin, dmax;
int openLeads;
long samples;
double normCnst = 32;
int current_out = 0;

void output(int out)
{
  sig4[current_out] = out;
  current_out++;
}

double detectQRS()
{
  
	double signal[BUFFSIZE], dcblock[BUFFSIZE], lowpass[BUFFSIZE], highpass[BUFFSIZE], derivative[BUFFSIZE], squared[BUFFSIZE], integral[BUFFSIZE], outputSignal[BUFFSIZE];
	int rr1[8], rr2[8], rravg1, rravg2, rrlow = 0, rrhigh = 0, rrmiss = 0;
	long unsigned int i, j, sample = 0, lastQRS = 0, lastSlope = 0, currentSlope = 0;
	int current;
	double peak_i = 0, peak_f = 0, threshold_i1 = 0, threshold_i2 = 0, threshold_f1 = 0, threshold_f2 = 0, spk_i = 0, spk_f = 0, npk_i = 0, npk_f = 0;
	bool qrs, regular = true, prevRegular;

	// Initializing the RR averages
	for (i = 0; i < 8; i++)
    {
        rr1[i] = 0;
        rr2[i] = 0;
    }

    do{
		if (sample >= BUFFSIZE)
		{
			for (i = 0; i < BUFFSIZE - 1; i++)
			{
				signal[i] = signal[i+1];
				dcblock[i] = dcblock[i+1];
				lowpass[i] = lowpass[i+1];
				highpass[i] = highpass[i+1];
				derivative[i] = derivative[i+1];
				squared[i] = squared[i+1];
				integral[i] = integral[i+1];
				outputSignal[i] = outputSignal[i+1];
			}
			current = BUFFSIZE - 1;
		}
		else
		{
			current = sample;
		}

        if (sample >= samples)
			break;

		signal[current] = sig0[sample];
		sample++; // Update sample counter

		// DC Block filter
		if (current >= 1)
			dcblock[current] = signal[current] - signal[current-1] + 0.995*dcblock[current-1];
		else
			dcblock[current] = 0;

		// Low Pass filter
		// y(nT) = 2y(nT - T) - y(nT - 2T) + x(nT) - 2x(nT - 6T) + x(nT - 12T)
		lowpass[current] = dcblock[current];
		if (current >= 1)
			lowpass[current] += 2 * lowpass[current-1];
		if (current >= 2)
			lowpass[current] -= lowpass[current-2];
		if (current >= 6)
			lowpass[current] -= 2 * dcblock[current-6];
		if (current >= 12)
			lowpass[current] += dcblock[current-12];

		// High Pass filter
		// y(nT) = 32x(nT - 16T) - [y(nT - T) + x(nT) - x(nT - 32T)]
		highpass[current] = -lowpass[current];
		if (current >= 1)
			highpass[current] -= highpass[current-1];
		if (current >= 16)
			highpass[current] += 32 * lowpass[current-16];
		if (current >= 32)
			highpass[current] += lowpass[current-32];

		// Derivative filter
		// y(nT) = (1/8T)[-x(nT - 2T) - 2x(nT - T) + 2x(nT + T) + x(nT + 2T)]
        derivative[current] = 0;
        if (current >= 2) // - x(nT - 2T) + x(nT + 2T)
        {
            derivative[current] -= 2 * highpass[current-2];
            derivative[current] += 2 * highpass[current+2];
        }
        if (current >= 1) // - 2x(nT - T) + 2x(nT + T)
        {
            derivative[current] -= 2 * highpass[current-1];
            derivative[current] += 2 * highpass[current+1];
            derivative[current] /= 8;
        }
        else
            derivative[current] = highpass[current];

		// Square function
		// y(nT) = [x(nT)]^2.
		squared[current] = derivative[current] * derivative[current];

		// Moving-Window Integration
		// y(nT) = (1/N)[x(nT - (N - 1)T) + x(nT - (N - 2)T) + ... x(nT)]

		integral[current] = 0;
		for (i = 0; i < WINDOWSIZE; i++)
		{
			if (current >= (double)i)
				integral[current] += squared[current - i];
			else
				break;
		}
		integral[current] /= (double)i;

		qrs = false;

		// If the current signal is above one of the thresholds (integral or filtered signal), it's a peak candidate
        if (integral[current] >= threshold_i1 || highpass[current] >= threshold_f1)
        {
            peak_i = integral[current];
            peak_f = highpass[current];
        }

		// If both the integral and the signal are above their thresholds, they're probably signal peaks
		if ((integral[current] >= threshold_i1) && (highpass[current] >= threshold_f1))
		{
			// There's a 200ms latency. If the new peak respects this condition, we can keep testing
			if (sample > lastQRS + FS / 5)
			{
			    // If it respects the 200ms latency, but it doesn't respect the 360ms latency, we check the slope
				if (sample <= lastQRS + (long unsigned int)(0.36 * FS))
				{
				    currentSlope = 0;
				    for (j = current - 10; j <= current; j++)
                        if (squared[j] > currentSlope)
                            currentSlope = squared[j];

				    if (currentSlope <= (double)(lastSlope / 2))
                    {
                        qrs = false;
                    }

                    else
                    {
                        spk_i = 0.125 * peak_i + 0.875 * spk_i;
                        threshold_i1 = npk_i + 0.25*(spk_i - npk_i);
                        threshold_i2 = 0.5*threshold_i1;

                        spk_f = 0.125 * peak_f + 0.875 * spk_f;
                        threshold_f1 = npk_f + 0.25 * (spk_f - npk_f);
                        threshold_f2 = 0.5 * threshold_f1;

                        lastSlope = currentSlope;
                        qrs = true;
                    }
				}
				// If it was above both thresholds and respects both latency periods, it's a R peak
				else
				{
				    currentSlope = 0;
                    for (j = current - 10; j <= current; j++)
                        if (squared[j] > currentSlope)
                            currentSlope = squared[j];

                    spk_i = 0.125 * peak_i + 0.875 * spk_i;
                    threshold_i1 = npk_i + 0.25 * (spk_i - npk_i);
                    threshold_i2 = 0.5 * threshold_i1;

                    spk_f = 0.125 * peak_f + 0.875 * spk_f;
                    threshold_f1 = npk_f + 0.25 * (spk_f - npk_f);
                    threshold_f2 = 0.5 * threshold_f1;

                    lastSlope = currentSlope;
                    qrs = true;
				}
			}
			// If the new peak doesn't respect the 200ms latency, it's noise. Update thresholds and move on to the next sample
			else
            {
                peak_i = integral[current];
				npk_i = 0.125 * peak_i + 0.875 * npk_i;
				threshold_i1 = npk_i + 0.25*(spk_i - npk_i);
				threshold_i2 = 0.5 * threshold_i1;
				peak_f = highpass[current];
				npk_f = 0.125 * peak_f + 0.875 * npk_f;
				threshold_f1 = npk_f + 0.25 * (spk_f - npk_f);
                threshold_f2 = 0.5*threshold_f1;
                qrs = false;
				outputSignal[current] = qrs;
				if (sample > BUFFSIZE)
                	output(outputSignal[0]);
                continue;
            }

		}

		// If a R-peak was detected, the RR-averages must be updated
		if (qrs)
		{
			// Add the newest RR-interval to the buffer and get the new average
			rravg1 = 0;
			for (i = 0; i < 7; i++)
			{
				rr1[i] = rr1[i+1];
				rravg1 += rr1[i];
			}
			rr1[7] = sample - lastQRS;
			lastQRS = sample;
			rravg1 += rr1[7];
			rravg1 *= 0.125;

			// If the newly-discovered RR-average is normal, add it to the "normal" buffer and get the new "normal" average
			// Update the "normal" beat parameters.
			if ( (rr1[7] >= rrlow) && (rr1[7] <= rrhigh) )
			{
				rravg2 = 0;
				for (i = 0; i < 7; i++)
				{
					rr2[i] = rr2[i+1];
					rravg2 += rr2[i];
				}
				rr2[7] = rr1[7];
				rravg2 += rr2[7];
				rravg2 *= 0.125;
				rrlow = 0.92 * rravg2;
				rrhigh = 1.16 * rravg2;
				rrmiss = 1.66 * rravg2;
			}

			prevRegular = regular;
			if (rravg1 == rravg2)
			{
				regular = true;
			}
			// If the beat had been normal but turned odd, change the thresholds
			else
			{
				regular = false;
				if (prevRegular)
				{
					threshold_i1 /= 2;
					threshold_f1 /= 2;
				}
			}
		}
		// If no R-peak was detected, it's important to check how long it's been since the last detection
		else
		{
		    // If no R-peak was detected for too long, use the lighter thresholds and do a back search
			if ((sample - lastQRS > (long unsigned int)rrmiss) && (sample > lastQRS + FS / 5))
			{
				for (i = current - (sample - lastQRS) + FS / 5; i < (long unsigned int)current; i++)
				{
					if ( (integral[i] > threshold_i2) && (highpass[i] > threshold_f2))
					{
					    currentSlope = 0;
                        for (j = i - 10; j <= i; j++)
                            if (squared[j] > currentSlope)
                                currentSlope = squared[j];

                        if ((currentSlope < (double)(lastSlope / 2)) && (i + sample) < lastQRS + 0.36 * lastQRS)
                        {
                            qrs = false;
                        }
                        else
                        {
                            peak_i = integral[i];
                            peak_f = highpass[i];
                            spk_i = 0.25 * peak_i+ 0.75 * spk_i;
                            spk_f = 0.25 * peak_f + 0.75 * spk_f;
                            threshold_i1 = npk_i + 0.25 * (spk_i - npk_i);
                            threshold_i2 = 0.5 * threshold_i1;
                            lastSlope = currentSlope;
                            threshold_f1 = npk_f + 0.25 * (spk_f - npk_f);
                            threshold_f2 = 0.5 * threshold_f1;
                            
                            //RR Average 1
                            rravg1 = 0;
                            for (j = 0; j < 7; j++)
                            {
                                rr1[j] = rr1[j+1];
                                rravg1 += rr1[j];
                            }
                            rr1[7] = sample - (current - i) - lastQRS;
                            qrs = true;
                            lastQRS = sample - (current - i);
                            rravg1 += rr1[7];
                            rravg1 *= 0.125;

                            //RR Average 2
                            if ( (rr1[7] >= rrlow) && (rr1[7] <= rrhigh) )
                            {
                                rravg2 = 0;
                                for (i = 0; i < 7; i++)
                                {
                                    rr2[i] = rr2[i+1];
                                    rravg2 += rr2[i];
                                }
                                rr2[7] = rr1[7];
                                rravg2 += rr2[7];
                                rravg2 *= 0.125;
                                rrlow = 0.92 * rravg2;
                                rrhigh = 1.16 * rravg2;
                                rrmiss = 1.66 * rravg2;
                            }

                            prevRegular = regular;
                            if (rravg1 == rravg2)
                            {
                                regular = true;
                            }
                            else
                            {
                                regular = false;
                                if (prevRegular)
                                {
                                    threshold_i1 /= 2;
                                    threshold_f1 /= 2;
                                }
                            }

                            break;
                        }
                    }
				}

				if (qrs)
                {
                    outputSignal[current] = false;
                    outputSignal[i] = true;
                    if (sample > BUFFSIZE)
                        output(outputSignal[0]);
                    continue;
                }
			}

			// Definitely no signal peak was detected.
			if (!qrs)
			{
				// If some kind of peak had been detected, then it's a noise peak. Thresholds must be updated accordinly.
				if ((integral[current] >= threshold_i1) || (highpass[current] >= threshold_f1))
				{
					peak_i = integral[current];
					npk_i = 0.125 * peak_i + 0.875 * npk_i;
					threshold_i1 = npk_i + 0.25 * (spk_i - npk_i);
					threshold_i2 = 0.5 * threshold_i1;
					peak_f = highpass[current];
					npk_f = 0.125 * peak_f + 0.875 * npk_f;
					threshold_f1 = npk_f + 0.25 * (spk_f - npk_f);
					threshold_f2 = 0.5*threshold_f1;
				}
			}
		}
        
		outputSignal[current] = qrs;
		if (sample > BUFFSIZE)
			output(outputSignal[0]);
	} while (sample < samples);

	// Output the last remaining samples on the buffer
	for (i = 1; i < BUFFSIZE; i++)
		output(outputSignal[i]);
}


int openRecord(char *record)
{
  return (isigopen(record, s, MAX_LEADS));
}

long ReadBuffer(long annot)
{
  long i,j;
  WFDB_Sample vec[MAX_LEADS];

  // A global variable s of type WFDB_Siginfo contains 
  // all signal information and is initialized using the 
  // isigopen (...) call (refer to function openRecord 
  // of this frame). 
  // Belos is an example of displaying the signal description, 
  // gain and baseline value for each signal in the opened 
  // signal file. 
  for (j=0;j<openLeads;j++)
    fprintf(stderr,"%s %lf %d\n", s[j].desc, s[j].gain, s[j].baseline);

  for (i=0;i<nrSamps;i++){
    if (getvec(vec)<0)
      break;
    sig0[i]=1000*(vec[0]-s[0].baseline)/s[0].gain;
    /* Transfer of signal into uV (this equals to multiplying signal 
       by 5), omit multiplication by 1000 to convert signal into mV. */    
  }
  return(nrSamps);
}

void writeQRS(char *record, int chan, char *ann)
{
  WFDB_Anninfo annIFO;
  WFDB_Annotation annot;
  int i;

  annIFO.name = ann; annIFO.stat = WFDB_WRITE;
  if (annopen(record, &annIFO, 1) < 0){
    fprintf(stderr,"Error opening QRS file\n");
    return;
  }
  annot.subtyp = annot.chan = annot.num = 0; annot.aux = NULL;
  for (i=0;i<samples;i++){
    if (sig4[i]!=0){
      annot.anntyp = sig4[i];
      annot.time=i;
      //printf("%d, %f\n", i, sig4[i]);
      if (putann(0, &annot) < 0) break;
    }
  }
}

void readQRS(char *record, char *ant, int chan)
{
  WFDB_Anninfo annIFO;
  WFDB_Annotation annot;
  long i;

  annIFO.name = ant; annIFO.stat = WFDB_READ;
  if (annopen(record, &annIFO, 1) < 0){
    fprintf(stderr,"Error opening QRS file\n");
    return;
  }

  annot.subtyp = annot.chan = annot.num = 0; annot.aux = NULL;
  for (i=0;i<samples;i++) 
    sig4[i]=0;

  i=0;
  while (getann(0, &annot)==0) 
  {
    if (annot.time>samples) 
    {
      fprintf(stderr,"Error reading annotation times\n");
      return;
    }
    sig4[annot.time]=annot.anntyp;
  }
}

int main(int argc, char *argv[])
{
  long i;
  char *record = NULL;
  char *annotator=NULL;
  
  for (i=1;i<argc;i++)
  {
    if (argv[i][0]!='-')
    {
      fprintf(stderr,"Error parsing command line\n");
      exit(1);
    }

    switch (argv[i][1]){
    case 'r':
      i++;
      record = (char *)malloc(sizeof(char)*(strlen(argv[i])+2));
      strcpy(record,argv[i]);
      break;

    case 'a': 
      i++;
      annotator = (char *) malloc(sizeof(char)*(strlen(argv[i])+2));
      strcpy(annotator, argv[i]);
      break;

    case 'n': 
      i++;
      normCnst = atof(argv[i]);
      break;

    default:
      fprintf(stderr,"Wrong switch\n");
      exit (2);
    }
  }

  if (record == NULL)
  {
    fprintf(stderr,"No record was specified, exiting\n");
    exit (2);
  }

  if ((openLeads=openRecord(record))<0)
  {
    fprintf(stderr,"Error opening record, exiting 1\n");
    exit (3);
  }

  fprintf(stderr,"Record opened\n");
  nrSamps = s->nsamp-1;
  sig0 = (double *)malloc(sizeof(double ) * nrSamps);
  sig4 = (double *)malloc(sizeof(double ) * nrSamps);

  if ((samples=ReadBuffer(0))<=0)
  {
    fprintf(stderr,"Error opening record, exiting\n");
    exit (3);
  }
  fprintf(stderr,"Data read\n");

  if (samples!=nrSamps)
    fprintf(stderr,"Sample count differs\n");

  detectQRS();
  fprintf(stderr,"Data analyzed\n");

  if (annotator==NULL)
  {
    annotator=(char *)malloc(sizeof(char)*5);
    sprintf(annotator,"qrs");
  }

  writeQRS(record, 4, annotator);
  fprintf(stderr,"Annotations written\n");
  wfdbquit();
  fprintf(stderr,"Record closed\n");

  free (sig0);
  free (sig4);
  free (annotator);

  return 0;
}
