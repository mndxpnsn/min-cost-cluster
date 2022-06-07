import java.util.Random;

public class Main {

    static double min(double x, double y) {
        double res = 0.0;

        if(x > 0) {
            if (x < y) res = x;
            else res = y;
        }

        return res;
    }

    static void init_mat(int n, double[][] mat) {
        // Initialize matrix mat with random data.
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                Random r = new Random();
                double min = 0.0;
                double max = 1.0;

                double rand_val = min + (max - min) * r.nextDouble();

                mat[i][j] = rand_val;
            }
        }
    }

    static double recur(double[][] X, int i, int j, double[][] dp) {
        double res = 0.0;

        if(dp[i][j] != 0.0) {
            return 0.0;
        }

        if(i == j - 1) {
            return X[i][j];
        }

        if(i < j - 1) {
            double min_val = 0.0;
            for(int k = i + 1; k < j - 1; ++k) {
                double val1 = recur(X, i, k, dp);
                double val2 = recur(X, k + 1, j, dp);
                min_val = min(min_val, val1 + val2);
            }

            res = X[i][j] + min_val;
        }

        dp[i][j] = res;

        return res;
    }

    static double  bottom_up(double[][] X, int n) {
        double res = 0.0;

        double[][] L = new double[n + 1][n + 1];

        for(int j = 1; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                if(i + j < n) {
                    double min_val = 0.0;
                    for(int k = i + 1; k < j - 1; ++k) {
                        min_val = min(min_val, L[i][k] + L[k + 1][j]);
                    }

                    L[i][i + j] = X[i][i + j] + min_val;
                }
            }
        }

        res = L[0][n - 1];

        return res;
    }

    static double top_down(double[][] X, int n) {
        double res = 0.0;

        double[][] dp = new double[n + 1][n + 1];

        res = recur(X, 0, n - 1, dp);

        return res;
    }

    public static void main(String[] args) {

        // Size of input vector
        int n = 5;

        // Declare and initialize cost matrix
        double[][] X = new double[n + 1][n + 1];

        init_mat(n, X);

        // Compute minimum cluster cost using a top-down strategy
        double cluster_cost1 = top_down(X, n);

        // Compute minimum cluster cost using a bottom-up strategy
        double cluster_cost2 = bottom_up(X, n);

        System.out.println("cluster_cost1: " + cluster_cost1);
        System.out.println("cluster_cost2: " + cluster_cost2);
    }
}
