pipeline {
    agent any
    stages {
        stage('Build') {
            steps {
                sh 'echo "I am working from"'
                sh '''
                    pwd
                '''
                sh '''
                    echo "Multiline shell steps works too"
                    ls -lah
                '''
            }
        }
    }
}