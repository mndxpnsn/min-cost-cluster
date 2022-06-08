import java.util.Random;

public class Main {

    static double min(double x, double y) {
        double res = 0.0;

        // Check for zeros
        if(x != 0.0 && y == 0.0) {
            return x;
        }

        if(x == 0.0 && y != 0.0) {
            return y;
        }

        // Compute minimum
        if (x < y) res = x;
        else res = y;

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

    static void print_mat(double[][] mat, int n) {
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                System.out.print(mat[i][j] + " ");
            }
            System.out.println("");
        }
    }

    static double recur(double[][] X, int i, int j, double[][] dp) {
        double res = 0.0;

        // Get results from memo table if available
        if(dp[i][j] != 0.0) {
            return dp[i][j];
        }

        // Return cluster of one node
        if(i == j - 1) {
            return X[i][j];
        }

        // Compute minimum cluster cost
        if(i < j - 1) {
            double min_val = 0.0;
            for(int k = i; k < j; ++k) {
                double val1 = recur(X, i, k, dp);
                double val2 = recur(X, k + 1, j, dp);
                min_val = min(val1 + val2, min_val);
            }

            res = X[i][j] + min_val;
        }

        // Store results in memo table
        dp[i][j] = res;

        return res;
    }

    static double bottom_up(double[][] X, int n) {
        double res = 0.0;

        double[][] L = new double[n + 1][n + 1];

        for(int j = 1; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                if(i + j < n) {
                    double min_val = 0.0;
                    for(int k = i; k <= i + j; ++k) {
                        min_val = min(L[i][k] + L[k + 1][i + j], min_val);
                    }

                    L[i][i + j] = X[i][i + j] + min_val;
                }
            }
        }

        res = L[0][n - 1];

        return res;
    }

    static double top_down(double[][] X, int n) {

        double[][] dp = new double[n + 1][n + 1];
        
        return recur(X, 0, n - 1, dp);
    }

    public static void main(String[] args) {

        // Size of input vector
        int n = 10;

        // Declare and initialize cost matrix
        double[][] X = new double[n + 1][n + 1];

        init_mat(n, X);

        // Compute minimum cluster cost using a top-down strategy
        double cluster_cost1 = top_down(X, n);

        // Compute minimum cluster cost using a bottom-up strategy
        double cluster_cost2 = bottom_up(X, n);

        System.out.println("cluster cost top-down: " + cluster_cost1);
        System.out.println("cluster cost bottom-up: " + cluster_cost2);
    }
}
