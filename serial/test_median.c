
int cmpfunc (double a, double b)
{
   return (a - b);
}

double sorted_median(double *vector, int n_items)
{
    qsort(vector, n_items, sizeof(double), cmpfunc);

    if(n_items % 2 != 0)
    {
        return vector[n_items/2];
    }
    else
    {
        return 0.5 * (vector[n_items/2] + [n_items/2 - 1]);
    }
}

double median(double *vector, int n_items)
{
    int full_splits = n_items/5;
    int semi_splits = n_items % 5;

    double *medians = (double *)malloc((full_splits + (semi_splits!=0 ? 1 : 0)) * sizeof(double));

    for(int i=0; i<full_splits; i++)
    {
        medians[i] = sorted_median(vector + 5*i, 5);
    }

    if(semi_splits != 0)
    {
        i++;
        medians[i] = sorted_median(vector + 5*i, semi_splits);
    }

    median(medians, full_splits + (semi_splits!=0 ? 1 : 0));

}