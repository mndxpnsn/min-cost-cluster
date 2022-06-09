import java.util.Random;

public class Main {

    static double LARGE_NUM = 1000.0;

    static double min(double x, double y) {
        double res = 0.0;

        // Check for zeros
        if(x != 0.0 && y == 0.0)
            return x;

        if(x == 0.0 && y != 0.0)
            return y;

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

    static void init_tree(int n, int[][] tree) {
        for(int i = 0; i < n - 1; ++i) {
            tree[i][i] = i;
            tree[i][i + 1] = i;
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

    static void print_tree_mat(int[][] tree, int n) {
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                System.out.print(tree[i][j] + " ");
            }
            System.out.println("");
        }
    }

    static double recur(double[][] mat, int i, int j, double[][] dp, int[][] tree) {
        double res = 0.0;

        // Get results from memo table if available
        if(dp[i][j] != 0.0) {
            return dp[i][j];
        }

        // Return cluster of one node
        if(i == j - 1) {
            tree[i][j] = i;
            return mat[i][j];
        }

        double bounds = LARGE_NUM;

        // Compute minimum cluster cost
        if(i < j - 1) {
            double min_val = 0.0;
            for(int k = i; k < j; ++k) {
                double val1 = recur(mat, i, k, dp, tree);
                double val2 = recur(mat, k + 1, j, dp, tree);
                min_val = min(val1 + val2, min_val);
                if(min_val < bounds) {
                    bounds = min_val;
                    tree[i][j] = k;
                }
            }

            res = mat[i][j] + min_val;
        }

        // Store results in memo table
        dp[i][j] = res;

        return res;
    }

    static double bottom_up(double[][] mat, int n, int[][] tree) {
        double res = 0.0;

        double[][] L = new double[n + 1][n + 1];

        for(int j = 1; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                if(i + j < n) {
                    double min_val = 0.0;
                    double bounds = LARGE_NUM;
                    for(int k = i; k <= i + j; ++k) {
                        min_val = min(L[i][k] + L[k + 1][i + j], min_val);
                        if(min_val < bounds) {
                            bounds = min_val;
                            tree[i][i + j] = k;
                        }
                    }

                    L[i][i + j] = mat[i][i + j] + min_val;
                }
            }
        }

        res = L[0][n - 1];

        return res;
    }

    static double top_down(double[][] mat, int n, int[][] tree) {

        double[][] dp = new double[n + 1][n + 1];
        
        return recur(mat, 0, n - 1, dp, tree);
    }

    static void print_tree_rec(int[][] tree, int i, int j, int n) {
        int k = tree[i][j];

        if(j - i == n - 1) {
            System.out.println("root node is: " + k + ", (" + i + "," + j + ")");
        }

        if(i == k && k + 1 < j) {
            int r = tree[k + 1][j];
            System.out.println("node " + r + " is the right child of node " + k + ", (" + (k + 1) + "," + j + ")");
            print_tree_rec(tree, k + 1, j, n);
        }
        
        if(k == j && i < k - 1) {
            int l = tree[i][k - 1];
            System.out.println("node " + l + " is the left child of node " + k + ", (" + i + "," + (k - 1) + ")");
            print_tree_rec(tree, i, k, n);
        }

        if(i < k) {
            int l = tree[i][k];
            System.out.println("node " + l + " is the left child of node " + k + ", (" + i + "," + k + ")");
            print_tree_rec(tree, i, k, n);
        }

        if(k + 1 < j && i < k) {
            int r = tree[k + 1][j];
            System.out.println("node " + r + " is the right child of node " + k + ", (" + (k + 1) + "," + j + ")");
            print_tree_rec(tree, k + 1, j, n);
        }
    }

    static void ver_tree_rec(double[][] mat, int[][] tree, int i, int j, int n, double[] ver_cost) {
        int k = tree[i][j];

        if(j - i == n - 1) {
            ver_cost[0] += mat[i][j];
        }

        if(i == k && k + 1 < j) {
            ver_cost[0] += mat[k + 1][j];
            ver_tree_rec(mat, tree, k + 1, j, n, ver_cost);
        }

        if(k == j && i < k - 1) {
            ver_cost[0] += mat[i][k - 1];
            ver_tree_rec(mat, tree, i, k, n, ver_cost);
        }

        if(i < k) {
            ver_cost[0] += mat[i][k];
            ver_tree_rec(mat, tree, i, k, n, ver_cost);
        }

        if(k + 1 < j && i < k) {
            ver_cost[0] += mat[k + 1][j];
            ver_tree_rec(mat, tree, k + 1, j, n, ver_cost);
        }
    }

    static void print_tree(int[][] tree, int n) {

        print_tree_rec(tree, 0, n - 1, n);
    }

    static double ver_tree(double[][] mat, int[][] tree, int n) {
        double[] ver_cost = new double[1];

        ver_tree_rec(mat, tree, 0, n - 1, n, ver_cost);

        return ver_cost[0];
    }

    public static void main(String[] args) {

        // Size of cluster tree + 1
        int n = 15;

        // Declare and initialize cost matrix and cluster trees
        double[][] cost = new double[n + 1][n + 1];
        int[][] ctree1 = new int[n + 1][n + 1];
        int[][] ctree2 = new int[n + 1][n + 1];

        init_mat(n, cost);

        // Compute minimum cluster cost using a top-down strategy
        double cluster_cost1 = top_down(cost, n, ctree1);

        // Compute minimum cluster cost using a bottom-up strategy
        double cluster_cost2 = bottom_up(cost, n, ctree2);

        // Verify construction of cluster tree obtained top-down
        double ver_cost1 = ver_tree(cost, ctree1, n);

        // Verify construction of cluster tree obtained bottom-up
        double ver_cost2 = ver_tree(cost, ctree2, n);

        System.out.println("tree top-down");
        print_tree_mat(ctree1, n);

        System.out.println("tree bottom-up");
        print_tree_mat(ctree2, n);

        System.out.println("cluster tree");
        print_tree(ctree1, n);

        System.out.println("cluster cost top-down: " + cluster_cost1);
        System.out.println("cluster cost bottom-up: " + cluster_cost2);
        System.out.println("Verification cluster cost top-down: " + ver_cost1);
        System.out.println("Verification cluster cost bottom-up: " + ver_cost2);
    }
}
